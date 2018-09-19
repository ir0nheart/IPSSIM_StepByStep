// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>
#include <iostream>
#include <thread>
#include <Windows.h>
#include <iomanip>
#include <vector>
#include <unordered_map>
#include <deque>
#include <fstream>
#include <string>
#include <streambuf>
#include <istream>
#include <viennacl/vector.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/gmres.hpp>
#include <viennacl/linalg/ilu.hpp>
#include <compcol_double.h>
#include <ilupre_double.h>
#include MATRIX_H
#include "gmres.h"
#include <ilupre_double.h>

#pragma comment(lib,"C:/Code Libraries/SparseLib++/1.7/lib/libsparse.lib")
#pragma comment(lib,"C:/Code Libraries/SparseLib++/1.7/lib/libspblas.lib")
#pragma comment(lib,"C:/Code Libraries/SparseLib++/1.7/lib/libmv.lib")
// TODO: reference additional headers your program requires here
