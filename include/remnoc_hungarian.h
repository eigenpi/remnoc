#ifndef _REMNOC_HUNGARIAN_H_
#define _REMNOC_HUNGARIAN_H_

#include <vector>


using namespace std;

////////////////////////////////////////////////////////////////////////////////
//
// HUNGARIAN_ONE
//
////////////////////////////////////////////////////////////////////////////////

// this is the (Minimum) Assignment Problem by Hungarian Algorithm; taken and
// adapted from Knuth's Stanford Graphbase; first size1 will be used for the
// number of IP/cores, and second size2 will be used for the candidate tiles
// - as a subset of tiles available on NOC -
// http://reptar.uta.edu/NOTES5311/hungarian.c
//
// initialize( dim1, dim2)
// set_cost( y, x) for any possible assignment
// run_hungarian() to solve the assignment problem
// report_assignment_of( y)

#define SIZE1_HUNGARIAN 64;
#define SIZE2_HUNGARIAN 64;
#define INFINITY_HUNGARIAN 0x7FFFFFFF;

class HUNGARIAN_ONE {
 private:
    int _size1;
    int _size2;
    vector<int> _array;
    vector<int> _result;
    // utils storage;
    vector<int> _unchosen_row; // _size1;
    vector<int> _row_mate; // _size2;
    vector<int> _row_dec;
    vector<int> _col_inc;
    vector<int> _parent_row;
    vector<int> _slack;
    vector<int> _slack_row;
 public:
    HUNGARIAN_ONE() {
        _size1 = SIZE1_HUNGARIAN;
        _size2 = SIZE2_HUNGARIAN;
        initialize( _size1, _size2);
    }
    ~HUNGARIAN_ONE() {}
    
    void initialize( int y, int x);
    void clear();
    void run_hungarian(); // solve it;
    int get_array( int i, int j) { return _array[i * _size2 + j]; }
    void set_array( int i, int j, int value) { _array[i * _size2 + j] = value; }
    void set_cost( int y, int x, int cost) { set_array( y, x, cost); }
    int report_assignment_of( int y) { return _result[y]; }
    void print_hungarian_assignment(); // for debug;
};

////////////////////////////////////////////////////////////////////////////////
//
// HUNGARIAN_TWO
//
////////////////////////////////////////////////////////////////////////////////


#endif
