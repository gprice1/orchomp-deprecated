import openravepy as r
import sys
import time
import numpy
import random


#globals
e = None
m = None 
robot = None
knownstates = None


def getRandomState():
    lower, upper = robot.GetActiveDOFLimits()
    value = []
    for i in xrange( len( upper ) ):
        value.append( lower[i] + (upper[i]-lower[i]) * random.random() )
    return value

def getRandomState():
    lower, upper = robot.GetDOFLimits()
    value = []
    for i in xrange( len( upper ) ):
        value.append( lower[i] + (upper[i]-lower[i]) * random.random() )
    return value
   
def getRandomStates( n ):
    return [ getRandomState() for i in xrange( n )]

def saveState( savefile, state ):
    savefile.write( "\t".join( [ str(i) for i in state] ) + "\n" )

def loadState( savefile ):
    line = savefile.readline()
    words = line.split()
    return [ float( i ) for i in words ]




def userGetRandomArmStates( filename ):


    getmore = True
    savefile = open( filename, 'w' )
    robot.SetActiveDOFs( robot.GetManipulators()[0].GetArmIndices() )
    
    for i in xrange( 100 ):
        
        state = []
        invalid = True
        while invalid:
            state = getRandomState( )
            robot.SetActiveDOFValues( state )
            
            for key, value in knownstates.items():
                invalid = False
                robot.SetDOFValues( value[1], robot.GetManipulators()[1].GetArmIndices() )
                
                if robot.CheckSelfCollision( ):
                    invalid = True
                    break
                
        state = [ state[i] for i in robot.GetManipulators()[0].GetArmIndices() ] 
        print state
        saveState( savefile, state )


def userGetRandomStates( filename ):

    getmore = True
    savefile = open( filename, 'w' )

    while getmore:

        while True:
            state = getRandomState()
            robot.SetDOFValues( state )
            if not robot.CheckSelfCollision( ):
                break

        while True:
            
            s = raw_input( '\n-->' )
            if s == 'q' or s == 'quit' or s == 'quit()' or s == 'end':
                getmore = False
                break

            elif s == 's' or s == 'save':
                saveState( savefile, state )
                break
            elif s =='n' or s=='no' or s=='next':
                break

            else:
                print "\nNot a valid command"

def readKnownStates( filename ):
    datafile = open( filename, 'r')
    data = datafile.read()
    data = data.replace( ',' , '')
    lines = data.split( "\n" )

    knownStates = {}

    for i in xrange( 0, len( lines ), 3 ):
        if i+2 >= len( lines ):
            break 
        tag = lines[i]
        tag = tag[:-1]

        state = []
        state.append( [ float( j ) for j in lines[i+1].split() ] )
        state.append( [ float( j ) for j in lines[i+2].split() ] )
        
        knownStates[ tag ] = state

    return knownStates


def setRightArmState( state ):
    robot.SetDOFValues( state, robot.GetManipulators()[0].GetArmIndices() )

def setstate( state ):

    robot.SetDOFValues( state[0], robot.GetManipulators()[0].GetArmIndices() )
    robot.SetDOFValues( state[1], robot.GetManipulators()[1].GetArmIndices() )

def main():
    global e
    global m
    global robot
    global knownstates

    e = r.Environment()
    m = r.RaveCreateModule( e, 'orcdchomp' )
    


    e.SetViewer( 'qtcoin' )
    #e.Load( 'lab1.env.xml' );

    e.Load( 'herb2_padded_nosensors.robot.xml')
    robot = e.GetRobot( 'Herb2' )
    
    knownstates = readKnownStates( "data/knownstates.data")
    userGetRandomArmStates( "data/randomarmstates_0.data" )

    return 
    
    
    #robot.SetActiveManipulator( 0 )
    
    filename = "randomstates2.data"
    #filename = "randomManipulatorStates_0.data"

    userGetRandomStates( filename )

    return

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
