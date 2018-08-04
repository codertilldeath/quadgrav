(ns quadgrav.quadtree)

(defrecord Point [x y])

(defprotocol IBoundary
  (inside? [this point])
  (overlap? [this other])
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
      (or (and (<= x ox         (+ x dx))
               (<= y oy         (+ y dy)))
          (and (<= x (+ ox odx) (+ x dx))
               (<= y (+ oy ody) (+ y dy)))
          (and (<= x ox         (+ x dx))
               (<= y (+ oy ody) (+ y dy)))
          (and (<= x (+ ox odx) (+ x dx))
               (<= y oy         (+ y dy)))))))

;; (every? #(overlap? (Boundary. 1 1 1 1) %1)
;;         [(Boundary. 1.5 1.5 1 1)
;;          (Boundary. 0.5 0.5 1 1)
;;          (Boundary. 1.5 0.5 1 1)
;;          (Boundary. 0.5 1.5 1 1)])

(defrecord Quadtree [bounds points thresh nw ne sw se])

(defn make-qtree [x y dx dy]
  (->Quadtree (->Boundary x y dx dy) [] 1 nil nil nil nil))

(defn subdivided? [qtree]
  (not (nil? (:quads qtree))))

(defn subdivide [qtree]
  (if (subdivided? qtree)
    qtree
    (let [{x :x
           y :y
           dx :dx
           dy :dy} (:bounds qtree)
          halfdx (/ dx 2)
          halfdy (/ dy 2)]
      (-> qtree
          (assoc :nw (make-qtree x            y            halfdx halfdy))
          (assoc :ne (make-qtree (+ x halfdx) y            halfdx halfdy))
          (assoc :sw (make-qtree x            (+ y halfdy) halfdx halfdy))
          (assoc :se (make-qtree (+ x halfdx) (+ y halfdy) halfdx halfdy))))))

(defn ins [qtree ref arr]
  (if (not (inside? (:bounds qtree) (:point @ref)))
    nil
    (let [n (count (:points qtree))]
      (if (< n (:thresh qtree))
        (conj arr :points)
        (or (ins (atom (:nw qtree)) ref (conj arr :nw))
            (ins (atom (:ne qtree)) ref (conj arr :ne))
            (ins (atom (:sw qtree)) ref (conj arr :sw))
            (ins (atom (:se qtree)) ref (conj arr :se)))))))

(defn insert [qtree ref]
  (update-in qtree (ins qtree ref []) #(conj % ref)))
  

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
    (insert (atom {:point (Point. 1.5 1.5)}))
    (insert (atom {:point (Point. 2 2)}))
    ;;(query (Boundary. 1 1 2 2 ))
    )
