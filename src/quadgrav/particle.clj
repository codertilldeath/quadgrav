(ns quadgrav.particle
  (:require [quadgrav.quadtree :refer [apply-forces orbital-velocity]]
            [quadgrav.point :as p]
            [quadgrav.point :refer [IPoint]])
  (:import [quadgrav.point Point])
)

(require '[taoensso.tufte :as tufte :refer (defnp p profiled profile)])

(defrecord Particle [mass point vel force on])

(defn make-particle [mass x y xvel yvel xforce yforce on]
  (->Particle mass
              (Point. x y)
              (if (= x 500)
                (Point. xvel yvel)
                (let [distance (p/sub (Point. x y) (Point. 500 500))
                      angle (+ (/ Math/PI 2) (.heading distance))]
                  (p/mul (Point. (Math/cos angle)
                                   (Math/sin angle))
                         (orbital-velocity 40000
                                           (+ (rand 10) (.mag distance))))))
              (Point. xforce yforce)
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
  (let [m (:mass particle)
        f (apply-forces tree theta particle)
        a (p/div f m)
        v (p/add a (:vel particle))
        p (p/add v (:point particle))]
    (-> particle
        (assoc :force f)
        (assoc :vel v)
        (assoc :point p))))


(defn make-particle-updater [tree theta xbound ybound]
  (fn [particle]
    (-> particle
;        (check-bounds xbound ybound)
        (update-forces tree theta))))


