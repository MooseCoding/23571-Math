# How to make a look up table for Decode

## Setup 

Make sure that you have the following installed with Pip:
- SciPY
- SciKit-Learn
- Matplotlib
- Numpy

```bash
pip install numpy scikit-learn matplotlib scipy
```

OR 

You can use the Google Colab link and copy it into your drive

[Colab](https://colab.research.google.com/drive/1qrTq7R3H0MQH1yPb6HKURfo0GvzGhuGH?usp=sharing), where the first cell is equivalent to [regression.py](https://github.com/MooseCoding/23571-Math/blob/main/regression.py) and the second cell is equivalent to [values.py](https://github.com/MooseCoding/23571-Math/blob/main/values.py). 

## 1. Manual Data Set

Find 6 common distances that you want to shoot at away from the center of the classifier (which is inside the classifier). My team calculated its about 6 inches inside the classifier from the location of the april tag of the classifier if you want to measure from that. 

Record the flywheel power and hood position that you used in order to obtain your desired accuracy (my team did 8 of 10 artifacts that made it in). Use different combinations of each from each distance. Note that you may have to have limitations on your flywheel power and/or your hood (e.g. a servo that can't extend beyond its 0.8 position). I know for a fact flywheel speed is constrained but servo position, if its an issue let me know on Discord and I'll rewrite it. 

## 2. Find Your Regression

Plug in your data in [regression.py](https://github.com/MooseCoding/23571-Math/blob/main/regression.py) where it says to put in your hood position, flywheel power and your constraints. Then run the program and the resultant graph would look great in a portfolio and the equation it spits out is important. You can change the polynomial degree but the builtin solver I provide in the other [values.py](https://github.com/MooseCoding/23571-Math/blob/main/values.py) only takes in a quadratic expression.

Your output should look like this: ![](https://github.com/MooseCoding/23571-Math/blob/main/IMG_0442.jpeg) 

## 3. Find Values For Distance

So plugin your values for a,b,c for the expression $ax^2 + bx + c = f(x)$, where x is the hood position into [values.py](https://github.com/MooseCoding/23571-Math/blob/main/values.py) Yes your b and c terms may involve y but that is fine. From there make sure that you plugin your start distance (recommended min is 6, which is as close to center as possible) and the end distance (recommended distance is 144 sqrt 2 since thats literally the maximum fire distance). Also specify the step size between the two distances. From there add in your flywheel constraints from earlier leave the (x,x,1000) part the same though. Then run the program and tada your lookup table. 

Your output should look like this: ![](https://github.com/MooseCoding/23571-Math/blob/main/IMG_0443.jpeg)

## 4. Find Home For Array

Copy and paste into any file of your choice in Java or Kotlin (and just use the auto switch feature). 

## 5. Table Lookup

From there you can write helper functions like I did to get an index of the lookup table from a distance. 

```kotlin
fun lerp(t:Double, low: DoubleArray, high: DoubleArray):DoubleArray {
    val hood = low[0] + t*(high[0] - low[0])
    val flywheel = low[1]+t*(high[1]-low[1])
    return doubleArrayOf(hood, flywheel)
}

fun getValues(distance:Double): DoubleArray {
    val high: Int = Math.ceil(distance).toInt()
    val low: Int = Math.floor(distance).toInt()
    return lerp(distance-low, points[low], points[high])
}
```

## 6. Implementation Ideas

Here's my implementation of this

```kotlin
  val dist: Double = sqrt((goalX-currentX).pow(2) + (goalY-currentY).pow(2)) // Current distance from the goal
  val values: DoubleArray = Aimbot.points[getIndex(dist)] // This is my lookup table

  Hood.updatePosition(values[0])  // The hood position
  Flywheels.updatePid(values[1]) // The target velocity 
```
