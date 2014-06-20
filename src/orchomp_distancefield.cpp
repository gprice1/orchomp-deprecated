#include "orchomp_distancefield.h"

namespace orchomp {


// a simple consturctor that initializes some values
DistanceField::DistanceField() : aabb_padding(0.2), cube_extent(0.02)
{
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
    std::string cube_name = "cube";
    OpenRAVE::KinBodyPtr cube = createCube( environment,
                                            pose_world,
                                            cube_name);

    //setup the grid object for computing the occupancy grid,
    //  and the gradient and distance fields
    setupGrid( sizes[0], sizes[1], sizes[2] );

    RAVELOG_INFO("computing occupancy grid ...\n");
    

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

        cube->SetTransform( world_to_cube );
        
        //check for collisions on the cube at the given transform
        if ( environment->CheckCollision( cube, kinbody )){
            //If there is a collision, set the value to 1
            grid(i,j,k) = 1;
            
            std::stringstream ss;
            ss << i << "_" << j << "_" << k;
            std::string nameofthecube = ss.str();
            createCube( environment, world_to_cube.trans, nameofthecube  );

        }else {
            //If there is no collision, set the value to -1.
            grid(i,j,k) = -1;
        }
    }
    }
    }
    
    //Creates a distance field and gradient field from the binary values.
    //  values above 0 are inside objects, 
    //  values below 0 are outside of objects
    grid.computeDistsFromBinary();
    
    //delete the cube , because we don't need this anymore
    environment->Remove( cube );
   
}

OpenRAVE::dReal DistanceField::getDist( const vec3 & trans, vec3 & gradient ){
    
    return grid.sample( trans, gradient );
}


} // namespace orchomp
