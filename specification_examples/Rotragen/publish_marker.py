from matplotlib.font_manager import json_dump
import rospy
from visualization_msgs.msg import Marker
import numpy as np
import json
f = open('robot.json')
jsondict = json.load(f)
rad = jsondict['global_params'][13]
obstacle = np.array(jsondict['global_params'])[range(10,13)]
rospy.init_node('rviz_marker')

marker_pub = rospy.Publisher("/visualization_marker", Marker, queue_size = 2)

marker = Marker()

marker.header.frame_id = "world"
marker.header.stamp = rospy.Time.now()

# set shape, Arrow: 0; Cube: 1 ; Sphere: 2 ; Cylinder: 3
marker.type = 2
marker.id = 0

# Set the scale of the marker
marker.scale.x = 2*rad 
marker.scale.y = 2*rad 
marker.scale.z = 2*rad 

# Set the color
marker.color.r = 0.0
marker.color.g = 1.0
marker.color.b = 0.0
marker.color.a = 1.0

# Set the pose of the marker
marker.pose.position.x = obstacle[0] 
marker.pose.position.y = obstacle[1] 
marker.pose.position.z = obstacle[2] 
marker.pose.orientation.x = 0.0
marker.pose.orientation.y = 0.0
marker.pose.orientation.z = 0.0
marker.pose.orientation.w = 1.0
print('radius ', rad)
print('origin ', obstacle)

marker_pub.publish(marker)
rospy.rostime.wallsleep(1.0)
marker_pub.publish(marker)