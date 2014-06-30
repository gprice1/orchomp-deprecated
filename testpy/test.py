import openravepy as r
import sys
import time
import numpy


local_commands = set([ 'environment' ])


def main():
    e = r.Environment()
    m = r.RaveCreateModule( e, 'orchomp' )
    e.SetViewer( 'qtcoin' )
    e.Load( 'lab1.env.xml' );

    # set the manipulator to leftarm
    #ikmodel = databases.inversekinematics.InverseKinematicsModel(
    #            robot,iktype=IkParameterization.Type.Transform6D)
    #if not ikmodel.load():
    #    ikmodel.autogenerate()

    Tz = r.matrixFromAxisAngle([0,0,numpy.pi/2])
    Tz[0,3] = 0.4  
    Tz[1,3] = 1.6
    
    print Tz
    with e:
        for body in e.GetBodies():
            body.SetTransform(numpy.dot(Tz,body.GetTransform()))


    time.sleep(3.0) #sleep for a while to allow the viewer to set up
    name = ''

    if len( sys.argv ) > 1:
        name = sys.argv[1]
    else:
        name = 'test_default.txt'
        
    f = open( name, 'r' )

    data = f.read()
    lines = data.split('\n')

    current_command = ''
    for line in lines : 
        
        '''for comments'''
        if len(line) == 0 or line[0] == '#':
            continue

        current_command += ' ' + line

        if line[-1] != "/":
            print "Command:" , current_command
            m.SendCommand( current_command )
            print "after command"
            current_command = ""
        else:
            current_command = current_command[:-1] + " "

   
    while True:
        s = raw_input( '-->' )
        if s == 'q' or s == 'quit' or s == 'quit()' or s == 'end':
            break
        m.SendCommand( s )

    del e

main()
