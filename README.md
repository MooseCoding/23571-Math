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

[Colab](https://colab.research.google.com/drive/1qrTq7R3H0MQH1yPb6HKURfo0GvzGhuGH?usp=sharing)

This is useful if you don't have administrator access to install python or any packages. 

## 1. Manual Data Set

Find 6-7 (or more) common distances that you want to shoot at away from the center of the classifier (which is inside the classifier). My team calculated its about 6 inches inside the classifier from the location of the april tag of the classifier if you want to measure from that. 

Record the flywheel power and hood position that you used in order to obtain your desired accuracy (my team did 8 of 10 artifacts that made it in). Use different combinations of each from each distance. Note that you may have to have limitations on your flywheel velocity and/or your hood (e.g. a servo that can't extend beyond its 0.8 position). These constraints are baked into the code already. Please make sure to accurately fill them out so you don't break something on your bot. 

Note: Please ensure you get at ONE sample point at the end points of your step size. E.g. if you want to tune touching the goal, please have a sample where you are touching the goal, same goes for far away. My recommended tuning spots are as follows: touching the classifier, then 4 points inside the close zone, 2 points in the far zone, and one point at the very far corner of the field (if you can shoot it). This is assuming you can shoot anywhere. 

## 2. Find Your Regression

Plug in your data in [tool.py](https://github.com/MooseCoding/23571-Math/blob/main/tool.py) where it says to put in your hood position, flywheel power and your constraints. Then run the program and the resultant graph would look great in a portfolio and the equation it spits out is important. It's what we are going to use to find our values for flywheel velocity and hood positioning. 

Your output should look like this: ![](https://github.com/MooseCoding/23571-Math/blob/main/IMG_0442.jpeg). 

## 3. Find Values For Distance

Baked into [tool.py](https://github.com/MooseCoding/23571-Math/blob/main/tool.py), we spit out your values for hood position and flywheel velocity for any given distance. It'll download a file to your computer called "Aimbot.java". This contains the array and all the helper functions you could possibly need. 

## 4. Find Home For Array

Shove this file anywhere in your repo. 

## 5. Table Lookup

Here are our built in functions that give you 

```kotlin
/**
*@param t: Double, a scale from 0-1 of how in between we are between points
*@param low: DoubleArray, the flywheel velocity and hood position at the closer point
*@param high: DoubleArray, the flywheel velocity and hood position at our further point
*@returns DoubleArray, which is the double array of flywheel velocity and hood position to shoot from at our current point
**/
fun lerp(t:Double, low: DoubleArray, high: DoubleArray):DoubleArray

/**
*@param distance: Double, the distance from the goal
*@returns DoubleArray, the flywheel velocity and hood position to shoot from
**/
fun getValues(distance:Double): DoubleArray
```

You will not directly interface with lerp unless you want to, I provide you with freedoms! Its called by default by `getValues()` but you can do whatever you want! Also: if the 
distance makes the array go out of bounds it will return the next closest values. Very useful for slipping odometry however I would always advise about having a manual control if things go awry with 
position tracking. 

## 6. Implementation Ideas

Here's my implementation of this

```kotlin
  val dist: Double = sqrt((goalX-currentX).pow(2) + (goalY-currentY).pow(2)) // Current distance from the goal
  val values: DoubleArray = Aimbot.getValues(dist)

  Hood.updatePosition(values[0])  // The hood position
  Flywheels.updatePid(values[1]) // The target velocity 
```

This is using NextFTC Subsystems but you can accomplish very similar things in any code base. 

## 7. Optional Forms

If you are using this tool this season please take a minute to fill out this survey for me to gather information on users to learn more (and also portfolio bait). It means a lot and even 
just filling this out does help me. So if you use the tool, I'd advise taking 2 minutes to put down your team number, years of coding experience, name, and an optional thank you message that will
be displayed in my portfolio <3. Link here: [Google Form](https://forms.gle/4BRSkeDXZFncZCJ47)

## 8. Extending this code

Awesumdude is going to use this for airsort adding/subtracting hood position and flywheel power depending on the sorting position. I advised him to do a linearly extrapolation of certain test values that result in the proper sorting technique for certain sample points. 

# Error Reporting

If you encounter any errors, please ping @MooseCoding in any FTC programming server (NextFTC, PSIKit, Dairy, #programming in FTC Discord). I will try to respond within the day, and I'll work through with you any issues you encounter!  
