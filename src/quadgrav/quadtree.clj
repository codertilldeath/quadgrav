(ns quadgrav.quadtree
  (:require [quil.core :as q]
            [taoensso.tufte :as tufte :refer (defnp p profiled profile)])
  (:import (processing.core PVector))
  )

(defrecord Point [x y])

(defn add [p1 p2]
  (let [{x :x
         y :y} p2]
    (-> p1
        (update :x (partial + x))
        (update :y (partial + y)))))

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
    (let [px (.x p)
          py (.y p)]
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
    (PVector. (+ x (/ dx 2))
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

(defn ins [qtree ref arr]
  (let [bounds (:bounds qtree)]
    (if (not (inside? bounds (:point @ref)))
      nil
      (let [n (count (:points qtree))]
        (if (< n (:thresh qtree))
          arr
          (let [qtree (subdivide qtree)]
            (or (ins (:nw qtree) ref (conj arr :nw))
                (ins (:ne qtree) ref (conj arr :ne))
                (ins (:sw qtree) ref (conj arr :sw))
                (ins (:se qtree) ref (conj arr :se)))))))))

(defn update-point-mass [oldpoint oldmass newpoint newmass]
  (/ (+ (* oldpoint oldmass)
        (* newpoint newmass))
     (+ oldmass newmass)))


(defn update-center [qtree ref]
  (let [mass (:mass qtree)
        newmass (:mass @ref)
        p (:center qtree)
        newp (:point @ref)]
    (set! (.x p) (update-point-mass (.x p) mass (.x newp) newmass))
    (set! (.y p) (update-point-mass (.y p) mass (.y newp) newmass))
    qtree))

(defn update-center-masses [qtree keys ref]
  (-> qtree
      (update-center ref)
      (update :mass #(+ % (:mass @ref)))
      (#(if-let [sub (first keys)]
          (update % sub (fn [x] (update-center-masses x (rest keys) ref)))
          %))))

(defn check-subdivide [qtree keys]
  (if (zero? (count keys))
    qtree
    (update-in qtree keys subdivide)))


(defnp insert [qtree ref]
  (let [keys (ins qtree ref [])]
    (-> qtree
        (check-subdivide keys)
        (update-in (conj keys :points) #(conj % ref))
        (update-center-masses keys ref))))

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

;; (defn query [qtree bounds]
;;   )



(overlap? (Boundary. 2.5 2 0.5 0.5) (Boundary. 1 1 2 2))

(def G (* 6.67408 (Math/pow 10 -1)))
;(def G (* 6.67408 (Math/pow 10 -11)))

(defn get-force [pv1 m1 pv2 m2]
  (let [grav (/ (* G m1 m2)
                (Math/pow (.dist pv1 pv2) 2))
        angle (.heading (.sub (.copy pv2) pv1))]
    (.mult (PVector. (Math/cos angle)
                     (Math/sin angle))
           grav)))

(defn barnes-hut-apply [tree theta sv pvector mass]
  (when tree
    (let [s (:dx (:bounds tree))
          d (.dist pvector (:center tree))]
      (when-not (zero? d)
        (if (> theta (/ s d))
          (.add sv (get-force pvector mass (:center tree) (:mass tree)))
          (do 
            (barnes-hut-apply (:nw tree) theta sv pvector mass)
            (barnes-hut-apply (:ne tree) theta sv pvector mass)
            (barnes-hut-apply (:sw tree) theta sv pvector mass)
            (barnes-hut-apply (:se tree) theta sv pvector mass)
            (doseq [point (:points tree)]
              (when-not (zero? (.dist pvector (:point @point)))
                (.add sv (get-force pvector mass (:point @point) (:mass @point)))))
            ))
        )
      sv)))

(defnp apply-forces [qtree theta particle]
  (let [pvector (:force @particle)]
    (.mult pvector 0)
    (barnes-hut-apply qtree theta pvector (:point @particle) (:mass @particle))))

(defn orbital-velocity [mass dist]
  (Math/sqrt (/ (* mass G) dist)))

(let [cljp1 (Point. 1 1)
      cljp2 (Point. 2 2)
      pv1 (PVector. 1 1)
      pv2 (PVector. 2 2)]
  (profile
   {}
   (dotimes [_ 1000]
     (p :pvector-add (.add pv1 pv2))
     (p :clojure-add (add cljp1 cljp2)))))
