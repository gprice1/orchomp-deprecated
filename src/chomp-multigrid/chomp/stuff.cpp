
class ArrayNode{
    size_t index;
    double value;
}

int array_size = num_spheres*2;
ArrayNode * nodes = new ArrayNode[ array_size ]; 

void assertSorted(){
    for ( int i = 0; i < array_size-1; i ++ ){
        assert( nodes[i].value < nodes[i+1].value );
    }
}

void sortArray(){
    for ( int i = 0; i < array_size-1; i ++ ){

        const ArrayNode current_node = nodes[i+1];

        //TODO should I remove this quick exit?
        if ( current_node.value > nodes[i].value ){ continue; }

        for ( int current_index = i;
              current_index >= 0;
              current_index-- )
        {
            if ( current_node.value < nodes[current_index].value ){
                nodes[current_index+1] = nodes[current_index];
            }
            else{
                nodes[current_index+1] = current_node;
                break;
            }
        }
    }
}

//this one is almost undoubtedly slower.
void sortArray2(){
    for ( int i = 0; i < array_size-1; i ++ ){

        const double current_value = nodes[i+1].value;
        if ( current_value > nodes[i].value ){ continue; }

        for ( int current_index = i;
              current_index >= 0;
              current_index-- )
        {
            //if the current element is in the correct position,
            //  push the other elements up, and insert the current 
            if ( current_value > nodes[current_index].value ){

                ArrayNode temp = nodes[i+1];
                
                //move everything over by one.
                for ( int j = i; j >= current_index; j --; ){
                    nodes[j+1] = nodes[j];
                }
                nodes[current_index] = temp;
                break;
            }
        }
    }
}


void getExtentsOfSDF( const DistanceField & df,
                      OpenRAVE::Vector lower,
                      OpenRAVE::Vector upper)
{
    
    lower = upper = df.world_to_grid.trans;
    
    OpenRAVE::Vector v = df.lengths;
    v = df.world_to_grid.rotate( v );
    for ( int i = 0; i < 3; i ++ ){
        if ( v[i] > 0 ){ upper[i] += v[i];}
        else{ lower[i] += v[i]; }
    }
}

