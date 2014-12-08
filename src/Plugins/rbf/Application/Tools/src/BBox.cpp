//-------------------------------------------------------------------
//
//  Permission is  hereby  granted, free  of charge, to any person
//  obtaining a copy of this software and associated documentation
//  files  ( the "Software" ),  to  deal in  the  Software without
//  restriction, including  without limitation the rights to  use,
//  copy, modify,  merge, publish, distribute, sublicense,  and/or
//  sell copies of the Software, and to permit persons to whom the
//  Software is  furnished  to do  so,  subject  to  the following
//  conditions:
//
//  The above  copyright notice  and  this permission notice shall
//  be included  in  all copies  or  substantial  portions  of the
//  Software.
//
//  THE SOFTWARE IS  PROVIDED  "AS IS",  WITHOUT  WARRANTY  OF ANY
//  KIND,  EXPRESS OR IMPLIED, INCLUDING  BUT NOT  LIMITED  TO THE
//  WARRANTIES   OF  MERCHANTABILITY,  FITNESS  FOR  A  PARTICULAR
//  PURPOSE AND NONINFRINGEMENT. IN NO EVENT  SHALL THE AUTHORS OR
//  COPYRIGHT HOLDERS  BE  LIABLE FOR  ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
//  ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
//  USE OR OTHER DEALINGS IN THE SOFTWARE.
//-------------------------------------------------------------------
//-------------------------------------------------------------------

#include <rbf/Application/Tools/src/BBox.h>

#include <limits>

BBox::BBox()
{
  empty = 1;
}

BBox::~BBox()
{
}

void BBox::setMin(vec3 _m)
{
  min = _m;
  isEmpty();
}

void BBox::setMin(float _m)
{
  min = vec3(_m,_m,_m);
  isEmpty();
}

void BBox::setMax(vec3 _m)
{
  max = _m;
  isEmpty();
}

void BBox::setMax(float _m)
{
  max = vec3(_m,_m,_m);
  isEmpty();
}

void BBox::reset()
{
  min = vec3(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
  max = vec3(std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min());
  isEmpty();
}

bool BBox::isEmpty()
{
  empty = !(min.x < max.x && min.y < max.y && min.z < max.z);
  return empty;
}

bool BBox::inside(const vec3 &pos)
{
  return (pos.x >= min.x && pos.x <= max.x &&
          pos.y >= min.y && pos.y <= max.y &&
          pos.z >= min.z && pos.z <= max.z);
}
