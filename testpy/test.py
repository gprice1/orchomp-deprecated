import openravepy as r
import sys
import time

def main():
    e = r.Environment()
    m = r.RaveCreateModule( e, 'orchomp' )
    e.SetViewer( 'qtcoin' )
    time.sleep(2.0) #sleep for a while to allow the viewer to set up
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
