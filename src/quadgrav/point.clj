(ns quadgrav.point)

(defprotocol IPoint
  (add [this other])
  (sub [this other])
  (mul [this scalar])
  (div [this scalar])
  (heading [this])
  (mag [this])
  (dist [this other]))

(defrecord Point [x y]
  IPoint
  (add [this other]
    (-> other
        (update :x #(+ x %))
        (update :y #(+ y %))))
  (sub [this other]
    (-> other
        (update :x #(- x %))
        (update :y #(- y %))))
  (mul [this scalar]
    (-> this
        (update :x #(* scalar %))
        (update :y #(* scalar %))))
  (div [this scalar]
    (-> this
        (update :x #(/ % scalar))
        (update :y #(/ % scalar))))
  (heading [this]
    (Math/atan2 y x))
  (mag [this]
    (let [sq #(* % %)]
      (Math/sqrt (+ (sq x)
                    (sq y)))))
  (dist [this other]
    (let [{ox :x
           oy :y} other
          sq #(* % %)]
      (Math/sqrt (+ (sq (- ox x))
                    (sq (- oy y)))))))
