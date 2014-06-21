
#include "orchomp_distancefield.h"

#define COLLISION -1
#define NOCOLLISION 1
#define XAXIS 0
#define YAXIS 1
#define ZAXIS 2


namespace orchomp {


// a simple consturctor that initializes some values
DistanceField::DistanceField() : aabb_padding(0.2), cube_extent(0.02)
{
}

//a simple test function to see if two grids are the same;
bool areEqual( DtGrid & first, DtGrid & second ){
    for ( int i = 0; i < first.nx(); i ++ ){
    for ( int j = 0; j < first.ny(); j ++ ){
    for ( int j = 0; k < first.nz(); k ++ ){
        if ( first(i,j,k) != second(i,j,k) ){
            return false;
        }
    }
    }
    }
}


//several collision routines
virtual inline bool DistanceField::isCollided( OpenRAVE::KinBodyPtr cube,
                                          const OpenRAVE::Transform & world_to_cube ){

    cube->SetTransform( world_to_cube );
    if ( environment->CheckCollision( cube, kinbody ) ){
        return true;
    }

    return false;
}

virtual bool DistanceField::isCollided( OpenRAVE::KinBodyPtr cube, int x1, int x2,
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

virtual bool DistanceField::isCollided( OpenRAVE::KinBodyPtr cube, int x, int y, int z )
{
    OpenRAVE::Vector offset( cube_extent + x*cube_extent*2.0, 
                             cube_extent + y*cube_extent*2.0,
                             cube_extent + z*cube_extent*2.0);

    OpenRAVE::Transform world_to_cube = pose_world * offset;

    return isCollided( cube, world_to_cube );

}
void DistanceField::setupGrid(size_t x, size_t y, size_t z )
{

    grid.clear();
    grid.resize( x, y, z, DtGrid::AXIS_Z,
                 cube_extent * 2,
                 vec3( pose_world[0],
                       pose_world[1],
                       pose_world[2] ));
}


//gets the transform from the origin to the center of the cell at
//  the specified indices.
inline void DistanceField::getCenterFromIndex( size_t x, size_t y, size_t z,
                                               OpenRAVE::Transform & t ) const {
    t.identity();
    const vec3 center = grid.cellCenter( x,y,z );
    t.trans[0] = center[0];
    t.trans[1] = center[1];
    t.trans[2] = center[2];
}

inline OpenRAVE::KinBodyPtr
DistanceField::createCube( OpenRAVE::EnvironmentBasePtr & env,
                           OpenRAVE::Vector & pos,
                           std::string & name)

{
    //create a cube to be used for collision detection in the world.
    //  create an object and name that object 'cube'
    OpenRAVE::KinBodyPtr cube = RaveCreateKinBody( env );
    cube->SetName( name.c_str() );
    //set the dimensions of the cube object
    std::vector< OpenRAVE::AABB > vaabbs(1);
    /* extents = half side lengths */
    vaabbs[0].extents = 
        OpenRAVE::Vector(cube_extent, cube_extent, cube_extent);
    vaabbs[0].pos = pos;
    cube->InitFromBoxes(vaabbs, 1);
    
    //add the cube to the environment
    env->Add( cube );

    return cube;

}


void DistanceField::createField( OpenRAVE::EnvironmentBasePtr & environment )
{
 

    //get the geometry information
    OpenRAVE::geometry::aabb< OpenRAVE::dReal > aabb;
    //OpenRAVE::KinBody::KinBodyStateSaver statesaver(kinbody);
    //kinbody->SetTransform(OpenRAVE::Transform());
    aabb = kinbody->ComputeAABB();
    
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
    
    //get the pose of the voxel grid with respect to the world frame:
    for ( size_t i = 0; i < 3; i ++ ){
        pose_world[i] = aabb.pos[i] - 0.5 * lengths[i];
    }

    RAVELOG_INFO("SDF pose: %f, %f , %f\n", pose_world[0],
                                            pose_world[1],
                                            pose_world[2]   );
    
    //create a cube to do collision detection
    std::string cube_name = "unitCube";
    unitcube = createCube( environment, pose_world, cube_name);


    //setup the grid object for computing the occupancy grid,
    //  and the gradient and distance fields
    setupGrid( sizes[0], sizes[1], sizes[2] );

    RAVELOG_INFO("computing occupancy grid ...\n");
    
    
    //fill the grid utilizing various methods.
    simplefill();
    octreefill( 0, grid.nx(), 0, grid.ny(), 0, grid.nz() );

    Bound b;
    b.bounds = { 0, grid.nx(), 0, grid.ny(), 0, grid.nz() };
    b.dists  = { grid.nx() - 1, grid.ny() - 1 , grid.nz() - 1};
    kdtreefill( b, XAXIS ); 
    //Creates a distance field and gradient field from the binary values.
    //  values above 0 are inside objects, 
    //  values below 0 are outside of objects
    grid.computeDistsFromBinary();
    
    //delete the cube , because we don't need this anymore
    environment->Remove( unitCube );
   
}

OpenRAVE::dReal DistanceField::getDist( const vec3 & trans, vec3 & gradient ){
    
    return grid.sample( trans, gradient );
}



OpenRAVE::KinBodyPtr DistanceField::createCube( int xdist, int ydist, int zdist )
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


void DistanceField::setgrid( int x1, int x2, int y1, int y2, int z1, int z2, int value ){
    for ( int i = x1; i < x2; i ++ ){
        for ( int j = y1; j < y2; j ++){
            for ( int k = z1; k < z2; k ++ ){

                grid( i, j, k ) = value;
            }
        }
    }
}


void DistanceField::simplefill(){

    //TODO make this do octree stuff
    //index over every point in the grid. and get the collision 
    //  information.
    for ( size_t i= 0; i < grid.nx() ; i ++ ){
    for ( size_t j= 0; j < grid.ny() ; j ++ ){
    for ( size_t k= 0; k < grid.nz() ; k ++ ){
        
        //find the transfrom from the grid to the center of the
        //  indexed cell
        OpenRAVE::Transform world_to_cube;
        getCenterFromIndex( i, j, k , world_to_cube);

        unitCube->SetTransform( world_to_cube );
        
        //check for collisions on the cube at the given transform
        if ( environment->CheckCollision( unitCube, kinbody )){
            //If there is a collision, set the value to 1
            grid(i,j,k) = COLLISION;
            
            //TODO remove this
            std::stringstream ss;
            ss << i << "_" << j << "_" << k;
            std::string nameofthecube = ss.str();
            createCube( environment, world_to_cube.trans, nameofthecube  );

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
    
    environment->Remove( cube )
}



void DistanceField::kdtreefill( Bound b, int AXIS ){

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

    environment->Remove( cube );
}

} // namespace orchomp
