int bounds[3];

#define COLLISION -1
#define NOCOLLISION 1
#define XAXIS 0
#define YAXIS 1
#define ZAXIS 2

typedef struct bounds_t {
    int bounds[6];
    int dists[3];

} Bound;


OpenRAVE::KinBodyPtr createCube( int xdist, int ydist, int zdist );
void setgrid( int x1, int x2, int y1, int y2, int z1, int z2, int value ){
inline isCollided( OpenRAVE::KinBodyPtr cube, const OpenRAVE::Transform & world_to_cube);


void kdtreefill( Bound b, int AXIS ){

    //Get the axis along which to split the space.
    for ( int i = 0; i < 3; i ++ ){
        if ( AXIS == ZAXIS ){ AXIS == XAXIS; }
        else { AXIS++ ; }

        //bail out when we find a dimension with greater than one dist.
        if (b.dist[AXIS] > 1 ){ break; }
    } 
    
    //base case : all of the Axes have a dist of 1, so it is just a single cube.
    if (b.dist[AXIS] == 1 )
    {
        if ( isCollided( unitCube, b.bounds[0], b.bounds[1], b.bounds[2] ) ){
            grid( b.bounds[0], b.bounds[0], b.bounds[0] ) = COLLISION;
        }else{
            grid( b.bounds[0], b.bounds[0], b.bounds[0] ) = NOCOLLISION;
        }
        return;
    }


    //If it is not a base case, split along the middle of the given dimension.
    b.dist[AXIS] = b.dist[AXIS]/2 + 1;
    const int midpoint = b.bounds[AXIS*2] + newdist;

    OpenRAVE::KinBodyPtr cube = createCube( b.dists[0], b.dists[1], b.dists[2] );

    /////////////Check the top half////////////////////////////////
    const int temp = b.bounds[AXIS*2+1];
    b.bounds[AXIS*2+1] = midpoint;

    if ( isCollided(cube, b.bounds[0], b.bounds[1], b.bounds[2],
                          b.bounds[3], b.bounds[4], b.bounds[5] ) )
    {
        kdtreefill( b, AXIS );
    }else {
        setgrid( b.bounds[0], b.bounds[1], b.bounds[2],
                 b.bounds[3], b.bounds[4], b.bounds[5], NOCOLLISION );
    }


    /////////////Check the bottom half////////////////////////////////
    //restore the original value of the upper bound;
    b.bounds[AXIS*2+1] = temp;
    b.bounds[AXIS*2] = midpoint;
    b.dist[AXIS] = b.bounds[AXIS*2+1] - b.bounds[AXIS*2];

    if ( isCollided( cube, b.bounds[0], b.bounds[1], b.bounds[2],
                           b.bounds[3], b.bounds[4], b.bounds[5] ) )
    {
        kdtreefill( b, AXIS );
    }else {
        setgrid( b.bounds[0], b.bounds[1], b.bounds[2],
                 b.bounds[3], b.bounds[4], b.bounds[5], NOCOLLISION );
    }

        
}


