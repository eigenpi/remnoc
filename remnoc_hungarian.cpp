#include "remnoc_hungarian.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//
// HUNGARIAN_ONE
//
////////////////////////////////////////////////////////////////////////////////

void HUNGARIAN_ONE::initialize( int y, int x) 
{
    _size1 = y;
    _size2 = x;
    clear(); // clear storage from previous use;
    
    int array_count = _size1 * _size2;
    int i = 0;
    _array.resize(array_count);
    for ( i = 0; i < array_count; i++) {
        _array[i] = -1;
    }
    assert( _array.size() == array_count);
    _result.resize( _size1);
    _unchosen_row.resize( _size1);
    _row_dec.resize( _size1);
    for ( i = 0; i < _size1; i++) {
        _result[i] = -1;
        _unchosen_row[i] = 0;
        _row_dec[i] = 0;
    }
    _row_mate.resize( _size2);
    _col_inc.resize( _size2);
    _parent_row.resize( _size2);
    _slack.resize( _size2);
    _slack_row.resize( _size2);
    for ( i = 0; i < _size2; i++) {
        _row_mate[i] = 0;
        _col_inc[i] = 0;
        _parent_row[i] = 0;
        _slack[i] = 0;
        _slack_row[i] = 0;
    }
}

void HUNGARIAN_ONE::clear() 
{
    _array.clear();
    _result.clear();
    _unchosen_row.clear();
    _row_mate.clear();
    _row_dec.clear();
    _col_inc.clear();   
    _parent_row.clear();
    _slack.clear();
    _slack_row.clear();
}

void HUNGARIAN_ONE::print_hungarian_assignment()
{
    // debug;
    cout << "\nTiles: ";
    for ( int j = 0; j < _size2; j++) {
        cout << "\t " << j;
    }
    for ( int i = 0; i < _size1; i++) {
        cout << "\nIP: " << i;
        for ( int j = 0; j < _size2; j++) {
            cout << "\t " << get_array(i, j);
        }
    }
    cout << "\nResult: ";
    for ( int i = 0; i < _size1; i++) {
        cout << "\nIP " << i << " --> " << _result[i];
    }
}

void HUNGARIAN_ONE::run_hungarian() 
{
    // solver;
    int j = 0, k = 0, l = 0, s = 0, t = 0, q = 0;
    int m = _size1;
    int n = _size2;
    int unmatched = 0;
    int cost = 0;

    // begin subtract column minima in order to start with lots of zeroes 12;
    for ( k = 0; k < m; k ++) {
        s = get_array(k, 0);
        for ( l = 1; l < n; l ++) {
            if ( get_array(k, l) < s) s = get_array(k, l);
        }
        if ( s) {
            for ( l = 0; l < n; l ++) {
                set_array(k, l, (get_array(k, l) - s));
            }
            cost += s;
        }
    }
    // begin initial state 16;
    t = 0;
    for ( l = 0; l < n; l ++) {
        _row_mate[l] = -1;
        _parent_row[l] = -1;
        _col_inc[l] = 0;
        _slack[l] = INFINITY_HUNGARIAN;
    }
    for ( k = 0; k < m; k ++) {
        s = 0;
        _row_dec[k] = s;
        for ( l = 0; l < n; l ++)
            if ((s == get_array(k, l)) && _row_mate[l] < 0) {
                _result[k] = l;
                _row_mate[l] = k;
                goto row_done;
            }
        _result[k] = -1;
        _unchosen_row[t++] = k;
      row_done:
        ;
    }
    
    // begin Hungarian algorithm 18;
    if ( !t)
        goto done;
    unmatched = t;
    while ( 1) {
        q = 0;
        while ( 1) {
            while ( q < t) {
                // begin explore node q of the forest 19;
                k = _unchosen_row[q];
                s = _row_dec[k];
                for ( l = 0; l < n; l ++) {
                    if ( _slack[l]) {
                        int del;
                        del = get_array(k, l) - s + _col_inc[l];
                        if ( del < _slack[l]) {
                            if ( !del) {
                                if (_row_mate[l] < 0) goto breakthru;
                                _slack[l] = 0;
                                _parent_row[l] = k;
                                _unchosen_row[t++] = _row_mate[l];
                            } else {
                                _slack[l] = del;
                                _slack_row[l] = k;
                            }
                        }
                    }
                }
                q++;
            }
            // begin introduce a new zero into the matrix 21;   
            s = INFINITY_HUNGARIAN;
            for ( l = 0; l < n; l ++) {
                if ( _slack[l] && _slack[l] < s) {
                    s = _slack[l];
                }
            }
            for ( q = 0; q < t; q ++) {
                _row_dec[_unchosen_row[q]] += s;
            }
            for ( l = 0; l < n; l ++) {
                if ( _slack[l]) {
                    _slack[l] -= s;
                    if ( !_slack[l]) {
                        // begin look at a new zero 22;
                        k = _slack_row[l];
                        if ( _row_mate[l] < 0) {
                            for ( j = l+1; j < n; j ++) {
                                if ( !_slack[j]) {
                                    _col_inc[j] += s; 
                                }
                            }
                            goto breakthru;
                        }
                        else {
                            _parent_row[l] = k;
                            _unchosen_row[t++] = _row_mate[l];
                        }
                    }
                } else {
                    _col_inc[l] += s;
                }
            }
        }
      breakthru:
        // begin update the matching 20;
        do {
            j = _result[k];
            _result[k] = l;
            _row_mate[l] = k;
            if (j < 0) break;
            k = _parent_row[j];
            l = j;
        } while (1);

        if ( --unmatched == 0) {
            goto done;
        }
        // begin get ready for another stage 17;
        t = 0;
        for ( l = 0; l < n; l ++) {
            _parent_row[l] =  -1;
            _slack[l] = INFINITY_HUNGARIAN;
        }
        for ( k = 0; k < m; k ++) {
            if ( _result[k] < 0) {
                _unchosen_row[t++] = k;
            }
        }
    }
  done:
    return;
}
