
#include "orchomp_distancefield.h"
#include "orchomp_mod.h"

#define COLLISION -2
#define NOCOLLISION -1
#define NOCOLLISION_EXPLORED 1

#define XAXIS 0
#define YAXIS 1
#define ZAXIS 2


namespace orchomp {


// a simple consturctor that initializes some values
DistanceField::DistanceField() : aabb_padding(0.2), cube_extent(0.02), 
    splitting_threshold( 100 )
{
}

//a simple test function to see if two grids are the same;
bool areEqual( DtGrid & first, DtGrid & second ){
    for ( size_t i = 0; i < first.nx(); i ++ ){
    for ( size_t j = 0; j < first.ny(); j ++ ){
    for ( size_t k = 0; k < first.nz(); k ++ ){
        if ( first(i,j,k) != second(i,j,k) ){
            return false;
        }
    }
    }
    }

    return true;
}


//several collision routines
inline bool DistanceField::isCollided( OpenRAVE::KinBodyPtr cube,
                              const OpenRAVE::Transform & world_to_cube ){

    timer.start( "collision" );
    cube->SetTransform( world_to_cube );

    bool collision = environment->CheckCollision( cube, kinbody );

    timer.stop( "collision" );

    return collision;
}

bool DistanceField::isCollided( OpenRAVE::KinBodyPtr cube, 
                                            int x1, int x2,
                                            int y1, int y2,
                                            int z1, int z2 )
{

    const int xdist = x2-x1;
    const int ydist = y2-y1;
    const int zdist = z2-z1;
    
    //get the offset from the grid to the center of the cube.
    OpenRAVE::Transform pose_grid_cube;
    
    pose_grid_cube.trans[0] = xdist*cube_extent + x1*cube_length; 
    pose_grid_cube.trans[1] = ydist*cube_extent + y1*cube_length;
    pose_grid_cube.trans[2] = zdist*cube_extent + z1*cube_length;

    OpenRAVE::Transform world_to_cube = pose_world_grid * pose_grid_cube;

    bool returnval =  isCollided( cube, world_to_cube );


    return returnval;
}

bool DistanceField::isCollided( OpenRAVE::KinBodyPtr cube,
                                            int x, int y, int z )
{
    
    

    OpenRAVE::Transform pose_grid_cube; 
    pose_grid_cube.trans[0] = cube_extent + x*cube_length;
    pose_grid_cube.trans[1] = cube_extent + y*cube_length;
    pose_grid_cube.trans[2] = cube_extent + z*cube_length;


    OpenRAVE::Transform world_to_cube = pose_world_grid * pose_grid_cube;
 
    bool retval = isCollided( cube, world_to_cube );


    return retval;
}


void DistanceField::setupGrid(size_t x, size_t y, size_t z )
{

    grid.clear();
    grid.resize( x, y, z, DtGrid::AXIS_Z,
                 cube_extent * 2,
                 vec3( 0,0,0));
}


//gets the transform from the origin to the center of the cell at
//  the specified indices.
inline void DistanceField::getCenterFromIndex( size_t x, size_t y, size_t z,
                                               OpenRAVE::Transform & t ) const {

    OpenRAVE::Transform pose_grid_cube; 
    pose_grid_cube.trans[0] = cube_extent + x*cube_length;
    pose_grid_cube.trans[1] = cube_extent + y*cube_length;
    pose_grid_cube.trans[2] = cube_extent + z*cube_length;
    
    //get the center of the cube in the world frame.
    //pose_world_center = pose_world_grid * pose_grid_center
    t = pose_world_grid * pose_grid_cube;
}

OpenRAVE::KinBodyPtr
DistanceField::createCube( OpenRAVE::EnvironmentBasePtr & env,
                           OpenRAVE::Transform & pos,
                           std::string & name)

{


    //create a cube to be used for collision detection in the world.
    //  create an object and name that object 'cube'
    OpenRAVE::KinBodyPtr cube = RaveCreateKinBody( env );

    if ( environment->GetKinBody( name.c_str() ).get() ){ return cube; }
    cube->SetName( name.c_str() );
    //set the dimensions of the cube object
    std::vector< OpenRAVE::AABB > vaabbs(1);
    /* extents = half side lengths */
    vaabbs[0].extents = 
        OpenRAVE::Vector(cube_extent, cube_extent, cube_extent);
    cube->InitFromBoxes(vaabbs, false );
    
    cube->SetTransform( pos );
    //add the cube to the environment
    env->Add( cube );

    return cube;

}


void DistanceField::createField( OpenRAVE::EnvironmentBasePtr & env )
{
    
    //set the cube length.
    cube_length = cube_extent * 2;
    
    //get the environment
    this->environment = env;

    OpenRAVE::geometry::aabb< OpenRAVE::dReal > aabb;
/* compute aabb when object is at world origin */
    {
        OpenRAVE::KinBody::KinBodyStateSaver statesaver(kinbody);
        kinbody->SetTransform(OpenRAVE::Transform());
        aabb = kinbody->ComputeAABB();
    }

    RAVELOG_INFO("    pos: %f %f %f\n", aabb.pos[0],
                                        aabb.pos[1],
                                        aabb.pos[2]);
    RAVELOG_INFO("extents: %f %f %f\n", aabb.extents[0],
                                        aabb.extents[1],
                                        aabb.extents[2]);
    
    size_t sizes[3];

    //compute the dimensions of the grid.
    //  size 0,1,2 are the number of cells
    //  in each x,y,z dimension.
    for ( int i = 0; i < 3; i ++ ){
        sizes[i] = size_t( ceil( (aabb.extents[i] + aabb_padding )
                                 / cube_extent ));
        
        RAVELOG_INFO("Grid Sizes [%d]: %d\n", i, sizes[i]);
    }
    
    //Compute the value of the lengths
    for ( int i = 0; i < 3; i ++ ){
        lengths[i] = sizes[i] * cube_extent * 2;
        RAVELOG_INFO("Grid Lengths [%d]: %f\n", i, lengths[i]);
    }
    
    OpenRAVE::Transform pose_kinbody_grid;

    //get the pose of the voxel grid with respect to the world frame:
    for ( size_t i = 0; i < 3; i ++ ){
        pose_kinbody_grid.trans[i] = aabb.pos[i] - 0.5 * lengths[i];
    }

    // pose_world_grid = pose_world_kinbody * pose_ kinbody_grid;
    pose_world_grid = kinbody->GetTransform() * pose_kinbody_grid;
    pose_grid_world = pose_world_grid.inverse();

    RAVELOG_INFO("SDF pose: %f, %f , %f\n", pose_world_grid.trans[0],
                                            pose_world_grid.trans[1],
                                            pose_world_grid.trans[2] );
    
    //create a cube to do collision detection
    std::string cube_name = "unitCube";
    unitCube = createCube( environment, pose_world_grid, cube_name);


    //setup the grid object for computing the occupancy grid,
    //  and the gradient and distance fields
    setupGrid( sizes[0], sizes[1], sizes[2] );

    size_t start, end[3];
    start = size_t( aabb_padding / cube_length );
    end[0] = grid.nx() - start;
    end[1] = grid.ny() - start;
    end[2] = grid.nz() - start;

    fillGridEdges( start, end );

    RAVELOG_INFO("computing occupancy grid ...\n");
    
    timer.start( "fill" );
    //fill the grid utilizing various methods.
    simplefill(start, end[0], start, end[1], start, end[2] );

    RAVELOG_INFO("Time to compute SimpleFill(): %f\n" , 
                 timer.stop( "fill" ) );

    RAVELOG_INFO("Time to compute collisions for SimpleFill(): %f\n" , 
                 timer.reset( "collision" ) );


    timer.start( "fill" );
    
    octreefill( start, end[0], start, end[1], start, end[2] );
    
    RAVELOG_INFO("Time to compute octreeFill(): %f\n" , 
                 timer.stop( "fill" ) );

    RAVELOG_INFO("Time to compute collisions for OctreeFill(): %f\n" , 
                 timer.reset( "collision" ) );

    
    timer.start( "fill" );
    
    Bound b(  start, end[0], start, end[1], start, end[2],
             end[0] - start, end[1] - start, end[2] - start );
    kdtreefill( b, XAXIS ); 
    
    RAVELOG_INFO("Time to compute KdtreeFill(): %f\n" , 
                 timer.stop( "fill" ) );

    RAVELOG_INFO("Time to compute collisions for KdtreeFill(): %f\n" , 
                 timer.reset( "collision" ) );


    
    simpleFloodFill( 0, grid.nx(), 0, grid.ny(), 0, grid.nz() );

    //delete the cube , because we don't need this anymore
    environment->Remove( unitCube );
   
}

OpenRAVE::dReal DistanceField::getDist( const OpenRAVE::Vector & pos,
                                        vec3 & gradient )
{

    //TODO remove some of this copying if you can.

    const OpenRAVE::Vector grid_point = pose_grid_world * pos;

    //check the bounds of the box, and if the point is not within the
    //  bounds, return huge_val;
    for ( size_t i = 0; i < 3; i ++){
        if ( grid_point[i] < 0 || grid_point[i] > lengths[i] ){
            return HUGE_VAL;
        }
    }

    vec3 trans( grid_point[0],
                grid_point[1],
                grid_point[2] );
    
    return grid.sample( trans, gradient );
}



OpenRAVE::KinBodyPtr DistanceField::createCube( int xdist, int ydist, int zdist )
{

    //create a cube to be used for collision detection in the world.
    //  create an object and name that object 'cube'
    OpenRAVE::KinBodyPtr cube = RaveCreateKinBody( environment );

    std::stringstream ss;
    ss   << kinbody->GetName() << "_"
         << xdist << "_"
         << ydist << "_"
         << zdist ; 

    const std::string name = ss.str();
    cube->SetName( name.c_str() );

    //set the dimensions of the cube object
    std::vector< OpenRAVE::AABB > vaabbs(1);

    /* extents = half side lengths */
    vaabbs[0].extents = 
        OpenRAVE::Vector(cube_extent*xdist,
                         cube_extent*ydist,
                         cube_extent*zdist );

    cube->InitFromBoxes(vaabbs, false);
    
    //add the cube to the environment
    environment->Add( cube );

    return cube;

}

//fills the reachable area with NOCOLLISION_EXPLORED
void DistanceField::simpleFloodFill(int x1, int x2, 
                                    int y1, int y2,
                                    int z1, int z2 ){

    vec3u lower( x1,y1,z1 ), upper( x2,y2,z2 );

    std::stack< vec3u > stack;
    stack.push( vec3u( x1,y1,z1 ) );

    assert( grid(x1,y1,z1) == NOCOLLISION );
    grid(x1,y1,z1) = NOCOLLISION_EXPLORED;

    while( !stack.empty() ){
        
        const vec3u top = stack.top();
        stack.pop();

        //for each of the three dimesions, try to grow out the viable area.
        for( int i = 0; i < 3; i ++ ){
            
            vec3u current = top;
            current[i] += 1;
            if ( current[i] < upper[i] && grid( current ) == NOCOLLISION ){
                grid( current ) = NOCOLLISION_EXPLORED;
                stack.push( current );
            }
            
            //since we already added 1, we need to minus 2 to get to
            //  top[i] - 1.
            current[i] -=2;
            if ( current[i] >= lower[i] && grid( current ) == NOCOLLISION ){
                grid( current ) = NOCOLLISION_EXPLORED;
                stack.push( current );
            }

        }
    }
        
}

void DistanceField::initVariableCube()
{

    //create a cube to be used for collision detection in the world.
    //  create an object and name that object 'cube'
    variableCube = RaveCreateKinBody( environment );


    const std::string name = "variableCube";
    variableCube->SetName( name.c_str() );

    //set the dimensions of the cube object
    variableCubeGeometry.resize(1);

    //add the cube to the environment
    environment->Add( variableCube );

}

void DistanceField::resizeVariableCube( int x, int y, int z ){

    variableCubeGeometry[0].extents[0] = x * cube_extent;
    variableCubeGeometry[0].extents[1] = y * cube_extent;
    variableCubeGeometry[0].extents[2] = z * cube_extent;

    variableCube->InitFromBoxes( variableCubeGeometry, false );

}

void DistanceField::setGrid( int x1, int x2, int y1, int y2, int z1, int z2, int value ){

    for ( int i = x1; i < x2; i ++ ){
        for ( int j = y1; j < y2; j ++){
            for ( int k = z1; k < z2; k ++ ){

                grid( i, j, k ) = value;
            }
        }
    }
}

//The padding creates a boundary around the object
//  that is not filled with anything, set all of these
//  areas to no collision.
//There is overlap in these areas, so do not double set them for
//  efficiency reasons
void DistanceField::fillGridEdges( size_t start, size_t * end ){
    
    //fill the z=0 face.
    setGrid( 0, grid.nx(), 0, grid.ny(), 0, start, NOCOLLISION );

    //fill the z=grid.nz() face.
    setGrid( 0, grid.nx(), 0, grid.ny(), end[2], grid.nz() , NOCOLLISION );

    //fill the y=0 face 
    setGrid( 0, grid.nx(), 0, start, start, end[2], NOCOLLISION);
    
    //fill the y =grid.ny face.
    setGrid( 0, grid.nx(), end[1], grid.ny(), start, end[2], NOCOLLISION );
    
    //fill the x = 0 face
    setGrid( 0, start, start, end[1], start, end[2], NOCOLLISION );

    //fill the x = grid.nx face
    setGrid( end[0], grid.nx(), start, end[1], start, end[2], NOCOLLISION );
}


void DistanceField::simplefill(size_t x1, size_t x2,
                               size_t y1, size_t y2,
                               size_t z1, size_t z2 ){
    
     
    for ( size_t i= x1; i < x2 ; i ++ ){
    for ( size_t j= y1; j < y2 ; j ++ ){
    for ( size_t k= z1; k < z2 ; k ++ ){
        
        //find the transfrom from the grid to the center of the
        //  indexed cell
        OpenRAVE::Transform world_to_cube;
        getCenterFromIndex( i, j, k , world_to_cube);

        timer.start( "collision" );
        unitCube->SetTransform( world_to_cube );
        //check for collisions on the cube at the given transform
        bool collision =  environment->CheckCollision( unitCube, kinbody );
        timer.stop( "collision" );


        if (collision ){
            //If there is a collision, set the value to 1
            grid(i,j,k) = COLLISION;
            
            /*
            //TODO remove this
            std::stringstream ss;
            ss << kinbody->GetName() << "_" << i << "_" << j << "_" << k;
            std::string nameofthecube = ss.str();
            createCube( environment, world_to_cube, nameofthecube  );
            */
        }else {
            //If there is no collision, set the value to -1.
            grid(i,j,k) = NOCOLLISION;
        }
    }
    }
    }
}

void DistanceField::octreefill( int x1, int x2, int y1, int y2, int z1, int z2 ){

    const int xdist = x2 - x1;
    const int ydist = y2 - y1;
    const int zdist = z2 - z1;
   
    if ( splitting_threshold > xdist * ydist * zdist )
    {
        simplefill( x1, x2, y1, y2, z1, z2);
        return;
    }

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
    
    const int new_xdist = (xdist-1)/2 + 1;
    const int new_ydist = (ydist-1)/2 + 1;
    const int new_zdist = (zdist-1)/2 + 1;
    
    const int xmid = x1 + new_xdist;
    const int ymid = y1 + new_ydist;
    const int zmid = z1 + new_zdist;
    

    //check all of the squares:

    OpenRAVE::KinBodyPtr cube = createCube( new_xdist,
                                            new_ydist,
                                            new_zdist);



    //check for collision in each of the 8 different spaces

    if( isCollided( cube, x1, xmid, y1, ymid, z1, zmid ) )
    {
        octreefill( x1, xmid, y1, ymid, z1, zmid );
    }else {
        setGrid( x1, xmid, y1, ymid, z1, zmid, NOCOLLISION );
    }

    if ( x_valid ){
        if( isCollided( cube, xmid, x2, y1, ymid, z1, zmid ) )
        {
            octreefill( xmid, x2, y1, ymid, z1, zmid );
        }else {
            setGrid( xmid, x2, y1, ymid, z1, zmid, NOCOLLISION );
        }
    }

    if ( y_valid ){
        if( isCollided( cube, x1, xmid, ymid, y2, z1, zmid ) )
        {
            octreefill( x1, xmid, ymid, y2, z1, zmid );
        }else {
            setGrid( x1, xmid, ymid, y2, z1, zmid, NOCOLLISION );
        }
    }

    if ( z_valid ){
        if( isCollided( cube, x1, xmid, y1, ymid, zmid, z2 ) )
        {
            octreefill( x1, xmid, y1, ymid, zmid, z2 );
        }else {
            setGrid( x1, xmid, y1, ymid, zmid, z2, NOCOLLISION );
        }
    }

    if ( x_valid && y_valid ){
        if( isCollided( cube, xmid, x2, ymid, y2, z1, zmid ) )
        {
            octreefill( xmid, x2, ymid, y2, z1, zmid );
        }else {
            setGrid( xmid, x2, ymid, y2, z1, zmid, NOCOLLISION );
        }
    }

    if ( y_valid && z_valid ){
        if( isCollided( cube, x1, xmid, ymid, y2, zmid, z2 ) )
        {
            octreefill( x1, xmid, ymid, y2, zmid, z2 );
        }else {
            setGrid( x1, xmid, ymid, y2, zmid, z2, NOCOLLISION );
        }
    }

    if ( z_valid && x_valid){
        if( isCollided( cube, xmid, x2, y1, ymid, zmid, z2 ) )
        {
            octreefill( xmid, x2, y1, ymid, zmid, z2 );
        }else {
            setGrid( xmid, x2, y1, ymid, zmid, z2, NOCOLLISION );
        }
    }

    if ( x_valid && y_valid && z_valid){
        if( isCollided( cube, xmid, x2, ymid, y2, zmid, z2 ) )
        {
            octreefill( xmid, x2, ymid, y2, zmid, z2 );
        }else {
            setGrid( xmid, x2, ymid, y2, zmid, z2, NOCOLLISION );
        }
    }
   


    environment->Remove( cube );

}



void DistanceField::kdtreefill( Bound b, int AXIS){
    
    if ( splitting_threshold > b.dists[0] * b.dists[1] * b.dists[2] )
    {
        simplefill( b.bounds[0], b.bounds[1],
                    b.bounds[2], b.bounds[3],
                    b.bounds[4], b.bounds[5] );
        return;
    }

    //Get the axis along which to split the space.
    for ( int i = 0; i < 3; i ++ ){
        if ( AXIS == ZAXIS ){ AXIS = XAXIS; }
        else { AXIS++ ; }

        //bail out when we find a dimension with greater than one dist.
        if (b.dists[AXIS] > 1 ){ break; }
    } 
    
    //base case : all of the Axes have a dist of 1,
    //  so it is just a single cube.
    if (b.dists[AXIS] == 1 )
    {
        if ( isCollided( unitCube, b.bounds[0], b.bounds[2], b.bounds[4] )){
            grid( b.bounds[0], b.bounds[2], b.bounds[4] ) = COLLISION;
        }else{
            grid( b.bounds[0], b.bounds[2], b.bounds[4] ) = NOCOLLISION;
        }
        return;
    }


    //If it is not a base case, split along the middle of the given
    // dimension.
    //  Round the dist up: this equation does that.
    //  ie: dist 4 -> 2, dist 3->2, dist 1->1
    b.dists[AXIS] = (b.dists[AXIS]-1)/2 + 1;
    const int midpoint = b.bounds[AXIS*2] + b.dists[AXIS];
    
    assert( midpoint < b.bounds[AXIS*2+1] );

    OpenRAVE::KinBodyPtr cube = createCube( b.dists[0],
                                            b.dists[1],
                                            b.dists[2] );




    /////////////Check the top half////////////////////////////////
    const int temp = b.bounds[AXIS*2+1];
    b.bounds[AXIS*2+1] = midpoint;

    if ( isCollided(cube, b.bounds[0], b.bounds[1], b.bounds[2],
                          b.bounds[3], b.bounds[4], b.bounds[5] ) )
    {
        kdtreefill( b, AXIS);
    }else {
        setGrid( b.bounds[0], b.bounds[1], b.bounds[2],
                 b.bounds[3], b.bounds[4], b.bounds[5], NOCOLLISION );
    }


    /////////////Check the bottom half////////////////////////////////
    //restore the original value of the upper bound;
    b.bounds[AXIS*2+1] = temp;
    b.bounds[AXIS*2] = midpoint;
    b.dists[AXIS] = b.bounds[AXIS*2+1] - b.bounds[AXIS*2];

    if ( isCollided( cube, b.bounds[0], b.bounds[1], b.bounds[2],
                           b.bounds[3], b.bounds[4], b.bounds[5] ) )
    {
        kdtreefill( b, AXIS);
    }else {
        setGrid( b.bounds[0], b.bounds[1], b.bounds[2],
                 b.bounds[3], b.bounds[4], b.bounds[5], NOCOLLISION );
    }

    environment->Remove( cube );

}

} // namespace orchomp
