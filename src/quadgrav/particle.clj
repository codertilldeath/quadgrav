(ns quadgrav.particle
  (:require [quadgrav.quadtree :refer [apply-forces orbital-velocity]])
  (:import (processing.core PVector)))


(require '[taoensso.tufte :as tufte :refer (defnp p profiled profile)])

(defrecord Particle [mass point vel force on])

(defn make-particle [mass x y xvel yvel xforce yforce on]
  (->Particle mass
              (PVector. x y)
              ;; (PVector. xvel yvel)
              (if (= x 500)
                (PVector. xvel yvel)
                (let [distance (.sub (PVector. x y) (PVector. 500 500))
                      angle (+ (/ Math/PI 2) (.heading distance))]
                  (.mult (PVector. (Math/cos angle)
                                   (Math/sin angle))
                         (orbital-velocity 20000
                                           (.mag distance)))))
              (PVector. xforce yforce)
              on))

(defnp check-bounds [particle xbound ybound]
  (let [pv (:point @particle)
        v (:vel @particle)]

    ;; Reflect the x axis accordingly
    (let [x (.x pv)
          vx (.x v)]
      (when (< x 0)
        (set! (.x v) (Math/abs vx)))
      
      (when (> x xbound)
        (set! (.x v) (- (Math/abs vx)))))

    
    ;; Reflect the y axis accordingly
    (let [y (.y pv)
          vy (.y v)]
      (when (< y 0)
        (set! (.y v) (Math/abs vy)))
      
      (when (> y ybound)
        (set! (.y v) (- (Math/abs vy))))))
  particle)

(defn update-forces [particle tree theta]
  (let [{f :force
         v :vel
         p :point
         m :mass} @particle]
    (apply-forces tree theta particle)
    (.div f m)
    (.add v f)
    (.add p v))
  particle)

(defn make-particle-updater [tree theta xbound ybound]
  (fn [particle]
    (-> particle
;        (check-bounds xbound ybound)
        (update-forces tree theta))))


