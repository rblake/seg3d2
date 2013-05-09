/**
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

// These are helper functions for the SDK samples (string parsing, timers, image helpers, etc)
#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

#ifdef WIN32
#pragma warning(disable:4996)
#endif

// includes, project
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <assert.h>
#include <math.h>

#include <fstream>
#include <vector>
#include <iostream>
#include <algorithm>

// includes, timer, string parsing, image helpers
#include <Application/BackscatterReconstruction/Algorithm/cuda_common/helper_timer.h>   // helper functions for timers
#include <Application/BackscatterReconstruction/Algorithm/cuda_common/helper_string.h>  // helper functions for string parsing
#include <Application/BackscatterReconstruction/Algorithm/cuda_common/helper_image.h>   // helper functions for image compare, dump, data comparisons
#include <Application/BackscatterReconstruction/Algorithm/cuda_common/exception.h>

#endif //  HELPER_FUNCTIONS_H
