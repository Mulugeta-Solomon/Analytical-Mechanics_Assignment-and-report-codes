robot = robotics.RigidBodyTree('DataFormat', 'column'); % robot configuration as column vector 

% add first body (fixed, with an offset from ICS)
body1 = robotics.RigidBody('body1');
body1.Joint = robotics.Joint('joint1', 'revolute');
T = trvec2tform([-0.5, 0 0.2])*eul2tform([pi/2 0 0]);
body1.Joint.setFixedTransform(T);
robot.addBody(body1, robot.BaseName);

% add second body (revolute joint)
body2 = robotics.RigidBody('body2');
body2.Joint = robotics.Joint('joint2', 'revolute');
body2.Joint.setFixedTransform([0.3 -pi/3 0.1 0], 'mdh');
robot.addBody(body2, 'body1');
