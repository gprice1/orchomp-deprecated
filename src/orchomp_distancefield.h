#ifndef _ORCHOMP_DISTANCEFIELD_H_
#define _ORCHOMP_DISTANCEFIELD_H_


#include "chomp-multigrid/mzcommon/DtGrid.h"
#include <openrave/openrave.h>
#include <openrave/planningutils.h>

namespace orchomp{

typedef vec3_t< OpenRAVE::dReal > vec3;


class DistanceField{

public:
    
    //the kinematic body that the Distance field branches from.
    OpenRAVE::KinBodyPtr kinbody;

    //aabb_padding : the padding for the bounding box,
    //cube_extent : half of the width/height/depth of a voxel cell
    //lengths : the height, width, and depth of the entire grid
    double aabb_padding, cube_extent, lengths[3];

    //the transform from the kinbody to the field.
    OpenRAVE::Transform pose_kinbody;

    //the transform from the world origin to the field
    OpenRAVE::Transform pose_world;
    
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
    //gets the transform from the origin to the center of the cell at
    //  the specified indices.
    void getCenterFromIndex( size_t x, size_t y, size_t z,
                           OpenRAVE::Transform & t ) const;    

    //sets up the DtGrid data structure to begin creating a 
    //  distance field
    void setupGrid(size_t x, size_t y, size_t z );
    
    inline OpenRAVE::KinBodyPtr createCube( OpenRAVE::EnvironmentBasePtr & env );

    void binaryFill();
    



};


} // orchomp namespace

#endif 
