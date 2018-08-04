(ns quadgrav.core
  (:require [quil.core :as q]
            [quil.middleware :as m]
            [quadgrav.quadtree :refer [make-qtree insert query qtree-size inside?]])
  (:import (quadgrav.quadtree Point Boundary)))

(defn update-particle [at]
  (swap! at #(update-in % [:point :x] (fn [x] (+ x (dec (rand 2))))))
  (swap! at #(update-in % [:point :y] (fn [x] (+ x (dec (rand 2))))))
  at)

(defrecord Particle [point on])

(defn new-state []
  (let [x (try
               (q/width)
               (catch Exception e 1200))
        y (try (q/height)
               (catch Exception e 800))]
    {:particles (repeatedly 10000 #(atom
                       (Particle. 
                        (Point. (rand x)
                                (rand y)
                         )
                        false)))
     :query (Boundary. 111 111 111 111)}))

(defn setup []
  (q/frame-rate 30)
  (q/color-mode :rgb)
  (new-state))

(defn update-state [state]
  (-> state
      (update-in [:particles] #(map update-particle %))))

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
  (let [{x :x
         y :y
         dx :dx
         dy :dy} (:bounds qtree)]
    (q/stroke 255)
    (q/no-fill)
    (q/stroke-weight 1)
    (q/rect x y dx dy))
  (if-let [quads (:quads qtree)]
    (doseq [[k v] quads]
      (draw-quadtree v))))

(defn draw-state [state]
  (try 
    (q/background 0)
    (let [particles (:particles state)
          tree (reduce #(insert %1 %2)
                       (make-qtree 0 0 (q/width) (q/height))
                       particles)
          inside (query tree (:query state))]
      (draw-quadtree tree)
      (q/stroke-weight 2)
      (q/stroke 0 255 0)
      (q/no-fill)
      (apply q/rect (vals (:query state)))
      (doseq [particle particles]
        (swap! particle #(update-in % [:on] (fn [x] false))))
      (doseq [particle inside]
        (swap! particle #(update-in % [:on] (fn [x] true))))
      (doseq [particle particles]
        (if (:on @particle)
          (q/fill 0 255 0)
          (q/fill 0 100 0))
        (q/stroke-weight 0)
        (let [pos (:point @particle)
              x (:x pos)
              y (:y pos)]
          ;;(q/with-translation [(/ (q/width) 2)
          ;;                     (/ (q/height) 2)]
          (q/ellipse x y 10 10))))
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

