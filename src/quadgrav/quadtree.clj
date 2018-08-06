(ns quadgrav.quadtree)

(defrecord Point [x y])

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
    (let [{px :x
           py :y} p]
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
    (let [bounds (:bounds qtree)]
      (-> qtree
          (assoc :nw (make-qtree (nw bounds)))
          (assoc :ne (make-qtree (ne bounds)))
          (assoc :sw (make-qtree (sw bounds)))
          (assoc :se (make-qtree (se bounds)))))))

(defn make-qtree
  ([x y dx dy]
   (subdivide (make-qtree (->Boundary x y dx dy))))
  ([bounds]
   (->Quadtree bounds [] 4 0 (center bounds) nil nil nil nil)))

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
        {newx :x
         newy :y} (:point @ref)]
    (-> qtree
        (update-in [:center :x] #(update-point-mass % mass newx newmass))
        (update-in [:center :y] #(update-point-mass % mass newy newmass)))))

(defn update-center-masses [qtree keys ref]
  (-> qtree
      (#(update-center % ref))
      (update :mass #(+ % (:mass @ref)))
      (#(if-let [sub (first keys)]
          (update % sub (fn [x] (update-center-masses x (rest keys) ref)))
          %))))


(defn insert [qtree ref]
  (let [keys (ins qtree ref [])]
    (-> qtree
        (#(if (zero? (count keys))
            %
            (update-in % keys subdivide)))
        (update-in (conj keys :points) #(conj % ref))
        (#(update-center-masses % keys ref)))))

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


(-> (make-qtree 1 1 2 2)
    (insert (atom {:point (Point. 1.5 1.5) :mass 100}))
    (insert (atom {:point (Point. 2 2) :mass 10}))
    (insert (atom {:point (Point. 2.1 2.1) :mass 10}))
    (insert (atom {:point (Point. 2.1 2.2) :mass 10}))
    (:center)
    )
(overlap? (Boundary. 2.5 2 0.5 0.5) (Boundary. 1 1 2 2))

(def G (* 6.67408 (Math/pow 10 -11)))

(defn query-force [qtree theta particle]
  {:x 0 :y (* theta 3)})
