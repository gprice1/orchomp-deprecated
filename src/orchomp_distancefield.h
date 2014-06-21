#ifndef _ORCHOMP_DISTANCEFIELD_H_
#define _ORCHOMP_DISTANCEFIELD_H_


#include "chomp-multigrid/mzcommon/DtGrid.h"
#include <openrave/openrave.h>
#include <openrave/planningutils.h>

namespace orchomp{

typedef vec3_t< OpenRAVE::dReal > vec3;


class DistanceField{

    //a small helper struct to keep track of bounds for the
    //  kd tree fill.
    typedef struct bounds_t {
        int bounds[6];
        int dists[3];

    } Bound;

public:
    
    //the kinematic body that the Distance field branches from.
    OpenRAVE::KinBodyPtr kinbody;

    OpenRAVE::KinBodyPtr unitCube;
    
    OpenRAVE::EnvironmentBasePrt environment;

    //aabb_padding : the padding for the bounding box,
    //cube_extent : half of the width/height/depth of a voxel cell
    //lengths : the height, width, and depth of the entire grid
    double aabb_padding, cube_extent, lengths[3];

    //the transform from the world origin to the field
    OpenRAVE::Vector pose_world;
    
    // The structure that holds the distance field, and the gradients
    typedef DtGrid_t<OpenRAVE::dReal> DtGrid;
    DtGrid grid;

    //PUBLIC FUNCTIONS:

    // a simple constructor that just initializes some values.
    DistanceField();
    
    //the main function that creates the distance field
    void createField( OpenRAVE::EnvironmentBasePtr & environment );
    
    ~DistanceField(){}

    OpenRAVE::dReal getDist( const vec3 & trans,
                             vec3 & gradient);
    

private:



    //various collision routines.
    virtual inline bool isCollided( OpenRAVE::KinBodyPtr cube,
                     const OpenRAVE::Transform & world_to_cube );
    virtual bool isCollided( OpenRAVE::KinBodyPtr cube,
                             int x1, int x2,
                             int y1, int y2,
                             int z1, int z2 );
    virtual bool DistanceField::isCollided( OpenRAVE::KinBodyPtr cube,
                                            int x, int y, int z );

    //gets the transform from the origin to the center of the cell at
    //  the specified indices.
    void getCenterFromIndex( size_t x, size_t y, size_t z,
                           OpenRAVE::Transform & t ) const;    

    //sets up the DtGrid data structure to begin creating a 
    //  distance field
    void setupGrid(size_t x, size_t y, size_t z );
    

    inline OpenRAVE::KinBodyPtr createCube(
                                OpenRAVE::EnvironmentBasePtr & env,
                                OpenRAVE::Vector & pos,
                                std::string & name);
    OpenRAVE::KinBodyPtr createCube( int xdist, int ydist, int zdist );

    void binaryFill();
    
    void setgrid( int x1, int x2, int y1, int y2,
                  int z1, int z2, int value );
    
    
    //fills the grid by simple iterating over every point and checking
    //  for collision.
    void simplefill();
    //two different methods of filling the binary collision grid.
    //Use a recursive octree method to split the space repeatedly into 8 
    //  squares
    void octreefill( int x1, int x2, int y1, int y2, int z1, int z2 );

    //use a recursive kdtree method to split space repeatedly in half. 
    void kdtreefill( Bound b, int AXIS );
};


} // orchomp namespace

#endif 
