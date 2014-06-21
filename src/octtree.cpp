
int bounds[3];

#define COLLISION -1
#define NOCOLLISION 1



OpenRAVE::KinBodyPtr createCube( int xdist, int ydist, int zdist )
{

    //create a cube to be used for collision detection in the world.
    //  create an object and name that object 'cube'
    OpenRAVE::KinBodyPtr cube = RaveCreateKinBody( env );

    std::stringstream ss;
    ss   << xdist << "_"
         << ydist << "_"
         << zdist ; 

    std::name = ss.str()
    cube->SetName( name.c_str() );

    //set the dimensions of the cube object
    std::vector< OpenRAVE::AABB > vaabbs(1);

    /* extents = half side lengths */
    vaabbs[0].extents = 
        OpenRAVE::Vector(cube_extent*xdist,
                         cube_extent*ydist,
                         cube_extent*zdist );

    cube->InitFromBoxes(vaabbs, 1);
    
    //add the cube to the environment
    env->Add( cube );

    return cube;

}


void setgrid( int x1, int x2, int y1, int y2, int z1, int z2, int value ){
    for ( int i = x1; i < x2; i ++ ){
        for ( int j = y1; j < y2; j ++){
            for ( int k = z1; k < z2; k ++ ){

                grid( i, j, k ) = value;
            }
        }
    }
}

inline isCollided( OpenRAVE::KinBodyPtr cube, const OpenRAVE::Transform & world_to_cube ){

    cube->SetTransform( world_to_cube );
    if ( environment->CheckCollision( cube, kinbody ) ){
        return true;
    }

    return false;
}

bool isCollided( OpenRAVE::KinBodyPtr cube, int x1, int x2,
                                            int y1, int y2,
                                            int z1, int z2 )
{
    const int xdist = x2-x1;
    const int ydist = y2-y1;
    const int zdist = z2-z1;
    
    //get the offset from the index to the center of the cube.
    OpenRAVE::Vector offset( xdist*cube_extent + x*cube_extent*2.0, 
                             ydist*cube_extent + y*cube_extent*2.0,
                             zdist*cube_extent + z*cube_extent*2.0);

    OpenRAVE::Transform world_to_cube = pose_world * offset;

    return isCollided( cube, world_to_cube );
}

bool isCollided( OpenRAVE::KinBodyPtr cube, int x, int y, int z )
{
    OpenRAVE::Vector offset( cube_extent + x*cube_extent*2.0, 
                             cube_extent + y*cube_extent*2.0,
                             cube_extent + z*cube_extent*2.0);

    OpenRAVE::Transform world_to_cube = pose_world * offset;

    return isCollided( cube, world_to_cube );

}

void recursiveOct( int x1, int x2, int y1, int y2, int z1, int z2 ){

    const int xdist = x2 - x1;
    const int ydist = y2 - y1;
    const int zdist = z2 - z1;
    
    //base case:
    if ( xdist == 1 && ydist == 1 && zdist == 1){
        if ( isCollided( unitCube, x1, y1, z1) ){
            grid( x1,y1,z1 ) = COLLISION;
        }else {
            grid( x1,y1,z1 ) = NOCOLLISION;
        }

        return;
    }
        
    //can we recurse over these dimensions?
    const bool x_valid = ( xdist > 1 ? true : false );
    const bool y_valid = ( ydist > 1 ? true : false );
    const bool z_valid = ( zdist > 1 ? true : false );
    
    const int xmid = x1 + xdist/2 + 1;
    const int xmid = y1 + ydist/2 + 1;
    const int xmid = z1 + zdist/2 + 1;
    

    //check all of the squares:
    cube = createcube( x1, xmid, y1, ymid, z1, zmid );

    
    //check for collision in each of the 8 different spaces

    if( isCollided( cube, x1, xmid, y1, ymid, z1, zmid ) )
    {
        recursiveOct( x1, xmid, y1, ymid, z1, zmid );
    }else {
        setGrid( x1, xmid, y1, ymid, z1, zmid, NOCOLLISION );
    }

    if ( x_valid ){
        if( isCollided( cube, xmid, x2, y1, ymid, z1, zmid ) )
        {
            recursiveOct( xmid, x2, y1, ymid, z1, zmid );
        }else {
            setGrid( xmid, x2, y1, ymid, z1, zmid, NOCOLLISION );
        }
    }

    if ( y_valid ){
        if( isCollided( cube, x1, xmid, ymid, y2, z1, zmid ) )
        {
            recursiveOct( x1, xmid, ymid, y2, z1, zmid );
        }else {
            setGrid( x1, xmid, ymid, y2, z1, zmid, NOCOLLISION );
        }
    }

    if ( z_valid ){
        if( isCollided( cube, x1, xmid, y1, ymid, zmid, z2 ) )
        {
            recursiveOct( x1, xmid, y1, ymid, zmid, z2 );
        }else {
            setGrid( x1, xmid, y1, ymid, zmid, z2, NOCOLLISION );
        }
    }

    if ( x_valid && y_valid ){
        if( isCollided( cube, xmid, x2, ymid, y2, z1, zmid ) )
        {
            recursiveOct( xmid, x2, ymid, y2, z1, zmid );
        }else {
            setGrid( xmid, x2, ymid, y2, z1, zmid, NOCOLLISION );
        }
    }

    if ( y_valid && z_valid ){
        if( isCollided( cube, x1, xmid, ymid, y2, zmid, z2 ) )
        {
            recursiveOct( x1, xmid, ymid, y2, zmid, z2 );
        }else {
            setGrid( x1, xmid, ymid, y2, zmid, z2, NOCOLLISION );
        }
    }

    if ( z_valid && x_valid){
        if( isCollided( cube, xmid, x2, y1, ymid, zmid, z2 ) )
        {
            recursiveOct( xmid, x2, y1, ymid, zmid, z2 );
        }else {
            setGrid( xmid, x2, y1, ymid, zmid, z2, NOCOLLISION );
        }
    }

    if ( x_valid && y_valid && z_valid){
        if( isCollided( cube, xmid, x2, ymid, y2, zmid, z2 ) )
        {
            recursiveOct( xmid, x2, ymid, y2, zmid, z2 );
        }else {
            setGrid( xmid, x2, ymid, y2, zmid, z2, NOCOLLISION );
        }
    }
    
    environment->Remove( cube )
}



