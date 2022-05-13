import numpy as np
import matplotlib.pyplot as plt


class AnglePlot:
    def __init__(self, line1, line2, offset=1.0, color=None, origin=(0, 0), len_x_axis=1, len_y_axis=1):
        self.line1 = line1
        self.line2 = line2
        self.offset = offset
        self.color = color
        self.origin = origin
        self.len_x_axis = len_x_axis
        self.len_y_axis = len_y_axis
    
    def angle_plot():
        xy1 = line1.get_xydata()
        xy2 = line2.get_xydata()
        slope1 = (xy1[1][1] - xy1[0][1]) / float(xy1[1][0] - xy1[0][0])
        angle1 = abs(math.degrees(math.atan(slope1)))
        slope2 = (xy2[1][1] - xy2[0][1]) / float(xy2[1][0] - xy2[0][0])
        angle2 = abs(math.degrees(math.atan(slope2)))
        theta1 = min(angle1, angle2)
        theta2 = max(angle1, angle2)
        self.angle = theta2 - theta1
        if color is None:
            color = line1.get_color()

        return patches.Arc(origin, len_x_axis * offset, len_y_axis * offset, 0, theta1, theta2, color=color, label=str(angle) + u"\u00b0")
    
    def get_angle():
        return self.angle
