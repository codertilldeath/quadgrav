(ns quadgrav.quadtree
  (:require [quil.core :as q]
            [taoensso.tufte :as tufte :refer (defnp p profiled profile)]
            [quadgrav.point :as p]
            [quadgrav.point :refer [IPoint]] )
  (:import quadgrav.point.Point)
  )

(defprotocol IBoundary
  (inside? [this point])
  (overlap? [this other])
  (nw [this])
  (ne [this])
  (sw [this])
  (se [this])
  (center [this])
  )

(defrecord Boundary [x y dx dy]
  IBoundary
  (inside? [this p]
    (let [px (:x p)
          py (:y p)]
      (and (<= x px)
           (< px (+ x dx))
           (<= y py)
           (< py (+ y dy)))))
  (overlap? [this o]
    (let [{ox :x
           oy :y
           odx :dx
           ody :dy} o]
      (not (or (< (+ ox odx) x)
               (> ox (+ x dx))
               (< (+ oy ody) y)
               (> oy (+ y dy))))))
  (nw [this]
    (->Boundary x              y              (/ dx 2) (/ dy 2)))
  (ne [this]
    (->Boundary (+ x (/ dx 2)) y              (/ dx 2) (/ dy 2)))
  (sw [this]
    (->Boundary x              (+ y (/ dy 2)) (/ dx 2) (/ dy 2)))
  (se [this]
    (->Boundary (+ x (/ dx 2)) (+ y (/ dy 2)) (/ dx 2) (/ dy 2)))
  (center [this]
    (Point. (+ x (/ dx 2))
            (+ y (/ dy 2)))))


;; (every? #(overlap? (Boundary. 1 1 1 1) %1)
;;         [(Boundary. 1.5 1.5 1 1)
;;          (Boundary. 0.5 0.5 1 1)
;;          (Boundary. 1.5 0.5 1 1)
;;          (Boundary. 0.5 1.5 1 1)])

(defrecord Quadtree [bounds points thresh mass center nw ne sw se])

(defn subdivided? [qtree]
  (not (nil? (:nw qtree))))

(declare make-qtree)

(defn subdivide [qtree]
  (if (subdivided? qtree)
    qtree
    (let [bounds (:bounds qtree)
          thresh (:thresh qtree)]
      (-> qtree
          (assoc :nw (make-qtree (nw bounds) thresh))
          (assoc :ne (make-qtree (ne bounds) thresh))
          (assoc :sw (make-qtree (sw bounds) thresh))
          (assoc :se (make-qtree (se bounds) thresh))))))

(defn make-qtree
  ([size x y dx dy]
   (subdivide (make-qtree (->Boundary x y dx dy) size)))
  ([x y dx dy]
   (subdivide (make-qtree (->Boundary x y dx dy) 4)))
  ([bounds size]
   (->Quadtree bounds [] size 0 (center bounds) nil nil nil nil)))

(defn update-point-mass [oldpoint oldmass newpoint newmass]
  (/ (+ (* oldpoint oldmass)
        (* newpoint newmass))
     (+ oldmass newmass)))

;; Change
(defn update-center [qtree ref]
  (let [mass (:mass qtree)
        newmass (:mass ref)
        newp (:point ref)]
    (-> qtree
        (update-in [:center :x]
                   #(update-point-mass % mass (:x newp) newmass))
        (update-in [:center :y]
                   #(update-point-mass % mass (:y newp) newmass)))))

(defn update-center-masses [qtree keys ref]
  (-> qtree
      (update-center ref)
      (update :mass #(+ % (:mass ref)))
      (#(if-let [sub (first keys)]
          (update % sub (fn [x] (update-center-masses x (rest keys) ref)))
          %))))

(defn check-subdivide [qtree keys]
  (if (zero? (count keys))
    qtree
    (update-in qtree keys subdivide)))

(defn ins [qtree ref arr]
  (let [bounds (:bounds qtree)]
    (if (not (inside? bounds (:point ref)))
      nil
      (let [n (count (:points qtree))]
        (if (< n (:thresh qtree))
          arr
          (let [qtree (subdivide qtree)]
            (or (ins (:nw qtree) ref (conj arr :nw))
                (ins (:ne qtree) ref (conj arr :ne))
                (ins (:sw qtree) ref (conj arr :sw))
                (ins (:se qtree) ref (conj arr :se)))))))))

(defnp insert [qtree ref]
  (let [keys (ins qtree ref [])]
    (-> qtree
        (check-subdivide keys)
        (update-in (conj keys :points) #(conj % ref))
        (update-center-masses keys ref))))

;; (-> (make-qtree 1 0 0 4 4)
;;     (insert {:point (Point. 0 0)
;;              :mass 10})
;;     (insert {:point (Point. 0 4)
;;              :mass 10})
;;     (insert {:point (Point. 4 0)
;;              :mass 20})
;;     (insert {:point (Point. 4 4)
;;              :mass 20})
;;     (:center))

(defn qtree-size [qtree]
  (if (nil? qtree)
    0
    (+ (count (:points qtree))
       (qtree-size (:nw qtree))
       (qtree-size (:ne qtree))
       (qtree-size (:sw qtree))
       (qtree-size (:se qtree)))))

(defn query [qtree bounds]
  (when qtree
    (let [tbounds (:bounds qtree)]
      (when (overlap? tbounds bounds)
        (concat (filter #(inside? bounds (:point @%))
                        (:points qtree))
                (query (:nw qtree) bounds)
                (query (:ne qtree) bounds)
                (query (:sw qtree) bounds)
                (query (:se qtree) bounds))))))


(overlap? (Boundary. 2.5 2 0.5 0.5) (Boundary. 1 1 2 2))

(def G (* 6.67408 (Math/pow 10 -1)))
;(def G (* 6.67408 (Math/pow 10 -11)))

(defn get-force [pv1 m1 pv2 m2]
  (let [grav (/ (* G m1 m2)
                (Math/pow (p/dist pv1 pv2) 2))
        angle (p/heading (p/sub pv2 pv1))]
    (p/mul (Point. (Math/cos angle)
                   (Math/sin angle))
         grav)))


(defn barnes-hut-apply [tree theta pvector mass]
  (when tree
    (let [s (:dx (:bounds tree))
          d (p/dist pvector (:center tree))]
      (when-not (zero? d)
        (if (> theta (/ s d))
          (get-force pvector mass (:center tree) (:mass tree))
          (do
            (reduce #(p/add %1 %2)
                    (Point. 0 0)
                    (filter #(not (nil? %))
                            (concat [(barnes-hut-apply (:nw tree) theta pvector mass)
                                     (barnes-hut-apply (:ne tree) theta pvector mass)
                                     (barnes-hut-apply (:sw tree) theta pvector mass)
                                     (barnes-hut-apply (:se tree) theta pvector mass)]
                                    (for [point (:points tree)]
                                      (if (< (p/dist pvector (:point point)) 10)
                                        (Point. 0 0)
                                        (get-force pvector mass (:point point) (:mass point)))))))))))))

(defnp apply-forces [qtree theta particle]
  (barnes-hut-apply qtree
                    theta
                    (:point particle)
                    (:mass particle)))

(defn orbital-velocity [mass dist]
  (Math/sqrt (/ (* mass G) dist)))

