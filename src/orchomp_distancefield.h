#ifndef _ORCHOMP_DISTANCEFIELD_H_
#define _ORCHOMP_DISTANCEFIELD_H_


#include "chomp-multigrid/mzcommon/DtGrid.h"
#include <openrave/openrave.h>
#include <openrave/planningutils.h>
#include "utils/timer.h"

namespace orchomp{

typedef vec3_t< OpenRAVE::dReal > vec3;

// The structure that holds the distance field, and the gradients
typedef DtGrid_t<OpenRAVE::dReal> DtGrid;

class DistanceField{

    //a small helper struct to keep track of bounds for the
    //  kd tree fill.
    typedef struct bounds_t {
        int bounds[6];
        int dists[3];
        bounds_t(int x1, int x2, int y1, int y2, int z1, int z2,
                 int xdist, int ydist, int zdist ){
            bounds[0] = x1;
            bounds[1] = x2;
            bounds[2] = y1;
            bounds[3] = y2;
            bounds[4] = z1;
            bounds[5] = z2;

            dists[0] = xdist;
            dists[1] = ydist;
            dists[2] = zdist;
        }
        std::string to_string(){
            std::stringstream ss;
            for ( int i = 0; i < 6; i ++ ){
                ss << bounds[i] << " ";
            }
            for ( int i = 0; i < 3; i ++ ){
                ss << dists[i] << " ";
            }
            std::string name = ss.str();
            return name;
        }
    } Bound;

public:
     
    //the kinematic body that the Distance field branches from.
    OpenRAVE::KinBodyPtr kinbody;

    OpenRAVE::KinBodyPtr unitCube, variableCube;
    std::vector< OpenRAVE::AABB > variableCubeGeometry;

    OpenRAVE::EnvironmentBasePtr environment;

    //aabb_padding : the padding for the bounding box,
    //cube_extent : half of the width/height/depth of a voxel cell
    //lengths : the height, width, and depth of the entire grid
    double aabb_padding, cube_extent, cube_length, lengths[3];

    //the transform from the world origin to the field
    OpenRAVE::Transform pose_world_grid, pose_grid_world;
   
    //the transform from the bottom left corner to the middle of the
    //cube
    OpenRAVE::Vector grid_center;

    void initVariableCube();
    void resizeVariableCube( int x, int y, int z );

    DtGrid grid;

    //PUBLIC FUNCTIONS:

    // a simple constructor that just initializes some values.
    DistanceField();
    
    //the main function that creates the distance field
    void createField( OpenRAVE::EnvironmentBasePtr & environment );
    
    ~DistanceField(){}

    OpenRAVE::dReal getDist( const OpenRAVE::Vector & pos ,
                             vec3 & gradient);
    
    int splitting_threshold;

private:

    Timer timer;

    //various collision routines.
    virtual inline bool isCollided( OpenRAVE::KinBodyPtr cube,
                     const OpenRAVE::Transform & world_to_cube );
    virtual bool isCollided( OpenRAVE::KinBodyPtr cube,
                             int x1, int x2,
                             int y1, int y2,
                             int z1, int z2 );
    virtual bool isCollided( OpenRAVE::KinBodyPtr cube,
                                            int x, int y, int z );

    //gets the transform from the origin to the center of the cell at
    //  the specified indices.
    void getCenterFromIndex( size_t x, size_t y, size_t z,
                           OpenRAVE::Transform & t ) const;    

    //sets up the DtGrid data structure to begin creating a 
    //  distance field
    void setupGrid(size_t x, size_t y, size_t z );
    

    OpenRAVE::KinBodyPtr createCube(
                                OpenRAVE::EnvironmentBasePtr & env,
                                OpenRAVE::Transform & pos,
                                std::string & name);

    OpenRAVE::KinBodyPtr createCube( int xdist, int ydist, int zdist );

    void binaryFill();
    
    void setGrid( int x1, int x2, int y1, int y2,
                  int z1, int z2, int value );
    
    void fillGridEdges( size_t start, size_t * end );
    
    //flood fill all of the reachable vaoxels in the grid.
    void simpleFloodFill( int x1, int x2, int y1, int y2, int z1, int z2 );

    //fills the grid by simple iterating over every point and checking
    //  for collision.
    void simplefill( size_t x1, size_t x2,
                     size_t y1, size_t y2,
                     size_t z1, size_t z2 );
    //two different methods of filling the binary collision grid.
    //Use a recursive octree method to split the space repeatedly into 8 
    //  squares
    void octreefill( int x1, int x2, int y1, int y2, int z1, int z2 );

    //use a recursive kdtree method to split space repeatedly in half. 
    void kdtreefill( Bound b, int AXIS);
};


} // orchomp namespace

#endif 
