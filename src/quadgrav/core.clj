(ns quadgrav.core
  (:require [quil.core :as q]
            [quil.middleware :as m]
            [quadgrav.quadtree :refer [make-qtree insert query qtree-size]]
            [quadgrav.particle :refer [make-particle make-particle-updater]])
  (:import (quadgrav.quadtree Boundary)
           (quadgrav.particle Particle)))

(require '[taoensso.tufte :as tufte :refer (defnp p profiled profile)])

(defnp new-state [particles qtree-depth theta]
  (let [x (try
               (q/width)
               (catch Exception e 1000))
        y (try (q/height)
               (catch Exception e 1000))
        halfx (/ x 2)
        halfy (/ y 2)
        quartx (/ halfx 2)
        quarty (/ halfy 2)
        particles (repeatedly particles #(atom
                       (make-particle 
                        (+ (rand 5) 5)
                        (+ quartx (rand halfx)) (+ quarty (rand halfy))
                        (dec (rand 2)) (dec (rand 2))
                        0 0
                        true)))]
    {:particles (conj particles (atom (make-particle 20000 500 500 0 0 0 0 true)))
     :query (Boundary. 111 111 111 111)
     :tree (reduce #(insert %1 %2)
                   (make-qtree 0 0 x y)
                   particles)
     :x x
     :y y
     :depth qtree-depth
     :theta theta}))

(defn setup []
  (q/frame-rate 30)
  (q/color-mode :rgb)
  (new-state 100 4 0.5))

(defn update-state [state]
  (-> state
      (assoc :tree (reduce #(insert %1 %2)
                           (make-qtree (:depth state) 0 0 (:x state) (:y state))
                           (:particles state)))
      (update-in [:particles] #(let [update-particle (make-particle-updater (:tree state)
                                                                            (:theta state)
                                                                            (:x state)
                                                                            (:y state))]
                                 (doall (map (fn [x] (update-particle x)) %))))))

(defn update-query [state event]
  (try
    (-> state
        (assoc :query
               (let [center-x (q/mouse-x)
                     center-y (q/mouse-y)]
                 (Boundary. (- center-x 50)
                            (- center-y 50)
                            100
                            100))))
    (catch Exception e)))

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
    (q/background 20)
    (let [particles (:particles state)
          tree (:tree state)
          inside (query tree (:query state))]
      ;;(draw-quadtree tree)
      (q/stroke-weight 2)
      (q/stroke 0 255 0)
      (q/no-fill)
      ;; (doseq [particle particles]
      ;;   (swap! particle #(update-in % [:on] (fn [x] false))))
      ;; (doseq [particle inside]
      ;;   (swap! particle #(update-in % [:on] (fn [x] true))))
      (apply q/rect (vals (:query state)))
      (doseq [particle particles]
        (if (:on @particle)
          (q/fill 255 255 255)
          (q/fill 0 100 0))
        (q/stroke-weight 0)
        (let [pos (:point @particle)
              x (int (.x pos))
              y (int (.y pos))
              mass (int (:mass @particle))]
          ;;(q/with-translation [(/ (q/width) 2)
          ;;                     (/ (q/height) 2)]
          (if (> mass 40)
            (q/ellipse x y 30 30)
            (q/ellipse x y mass mass)))))
    (catch Exception e)))

(defn -main [& args]
  (q/defsketch quadgrav
    :title "You spin my circle right round"
    :size [1000 1000]
    :setup setup
    :update update-state
    :mouse-moved update-query
    :draw draw-state
    :features [:keep-on-top]
    :middleware [m/fun-mode]))


(tufte/add-basic-println-handler! {})

(let [a (new-state 100 20 1)]
  (profile
   {}
   (let [b (update-state a)])))
   

