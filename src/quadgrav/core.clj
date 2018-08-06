(ns quadgrav.core
  (:require [quil.core :as q]
            [quil.middleware :as m]
            [quadgrav.quadtree :refer [make-qtree insert query qtree-size query-force]])
  (:import (quadgrav.quadtree Point Boundary)))

(defn check-bounds [at]
  (let [{x :x
         y :y} (:point @at)]
    (when (< x 0)
      (swap! at #(update-in % [:vel :x] (fn [x] (Math/abs x)))))
    (when (> x (q/width))
      (swap! at #(update-in % [:vel :x] (fn [x] (- (Math/abs x))))))
    (when (< y 0)
      (swap! at #(update-in % [:vel :y] (fn [x] (Math/abs x)))))
    (when (> y (q/height))
      (swap! at #(update-in % [:vel :y] (fn [x] (- (Math/abs x))))))))

(defn get-force [part tree]
  (let [{x :x
         y :y} (query-force tree 0.5 part)]
    (Point. x y)))

(defn update-particle [at tree]
  (check-bounds at)
  (swap! at #(assoc % :force (get-force % tree)))
  (swap! at #(update-in % [:vel :x] (fn [x] (+ x (/ (:x (:force %)) (:mass %))))))
  (swap! at #(update-in % [:vel :y] (fn [x] (+ x (/ (:y (:force %)) (:mass %))))))
  (swap! at #(update-in % [:point :x] (fn [x] (+ x (:x (:vel %))))))
  (swap! at #(update-in % [:point :y] (fn [x] (+ x (:y (:vel %))))))
  at)

(defrecord Particle [mass point vel force on])

(defn new-state []
  (let [x (try
               (q/width)
               (catch Exception e 1200))
        y (try (q/height)
               (catch Exception e 800))
        particles (repeatedly 200 #(atom
                       (Particle. 
                        (+ (rand 5) 5)
                        (Point. (rand x)
                                (rand y))
                        (Point. (dec (rand 2))
                                (dec (rand 2)))
                        (Point. 0 4)
                        false)))]
    {:particles particles
     :query (Boundary. 111 111 111 111)
     :tree (reduce #(insert %1 %2)
                   (make-qtree 0 0 (q/width) (q/height))
                   particles)}))

(defn setup []
  (q/frame-rate 30)
  (q/color-mode :rgb)
  (new-state))

(defn update-state [state]
  (-> state
      (assoc :tree (reduce #(insert %1 %2)
                           (make-qtree 0 0 (q/width) (q/height))
                           (:particles state)))
      (update-in [:particles] #(map (fn [x] (update-particle x (:tree state))) %))))

(defn update-query [state event]
  (-> state
      (assoc :query
             (let [center-x (q/mouse-x)
                   center-y (q/mouse-y)]
               (Boundary. (- center-x 50)
                          (- center-y 50)
                          100
                          100)))))

(defn draw-quadtree [qtree]
  (when qtree
    (let [{x :x
           y :y
           dx :dx
           dy :dy} (:bounds qtree)]
      (q/stroke 255)
      (q/no-fill)
      (q/stroke-weight 1)
      (q/rect x y dx dy))
    (draw-quadtree (:nw qtree))
    (draw-quadtree (:ne qtree))
    (draw-quadtree (:sw qtree))
    (draw-quadtree (:se qtree))))

(defn draw-state [state]
  (try 
    (q/background 0)
    (let [particles (:particles state)
          tree (:tree state)
          inside (query tree (:query state))]
      (draw-quadtree tree)
      (q/stroke-weight 2)
      (q/stroke 0 255 0)
      (q/no-fill)
      (doseq [particle particles]
        (swap! particle #(update-in % [:on] (fn [x] false))))
      (doseq [particle inside]
        (swap! particle #(update-in % [:on] (fn [x] true))))
      (apply q/rect (vals (:query state)))
      (doseq [particle particles]
        (if (:on @particle)
          (q/fill 0 255 0)
          (q/fill 0 100 0))
        (q/stroke-weight 0)
        (let [pos (:point @particle)
              x (:x pos)
              y (:y pos)
              mass (:mass @particle)]
          ;;(q/with-translation [(/ (q/width) 2)
          ;;                     (/ (q/height) 2)]
          (q/ellipse x y mass mass))))
    (catch Exception e)))

(defn -main [& args]
  (q/defsketch quadgrav
    :title "You spin my circle right round"
    :size [1200 800]
    :setup setup
    :update update-state
    :mouse-moved update-query
    :draw draw-state
    :features [:keep-on-top]
    :middleware [m/fun-mode]))

