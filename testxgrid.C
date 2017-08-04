#include "XGridUtil.h"
using namespace ESMCI;
int main(int argc, char * argv[]){
            { // One edge of the cut polygon overlaps with the subject polygon internally on the right
              double p[8] = {58, -19.8, 56, -19.8, 56, -18, 58, -18 };
              double q[8] = {57.85, -18.65, 57.5, -18.65, 57.5, -18.2715, 57.85, -18.2715};
              test_clip2D(2,3,4,p,4,q);
              std::vector<polygon> differences;
              weiler_clip_difference(2,2,4,p,4,q, differences);
            } 
  if(0){
            { // s contains c with one common vertex
              double p[6] = {72.043,  18.1369, 73.9063,  16.6181, 72.5,  16.8274};
              double q[8] = {73,  16.8274, 72.5,  16.8274, 72.5,  17.178, 73,  17.178};
              test_clip2D(2,3,3,p,4,q);
            } 
            { // One edge of the cut polygon overlaps with the subject polygon internally on the left
              double p[8] = {36, -18, 34, -18, 34, -16.2, 36, -16.2};
              double q[8] = {34.5, -16.8274, 34, -16.8274, 34.2, -16.4836, 34.5, -16.4836};
              test_clip2D(2,3,4,p,4,q);
              std::vector<polygon> differences;
              weiler_clip_difference(2,2,4,p,4,q, differences);
            } 
            { // One edge of the cut polygon overlaps with the subject polygon internally on the right
              double p[8] = {58, -19.8, 56, -19.8, 56, -18, 58, -18 };
              double q[8] = {58, -18.65, 57.5, -18.65, 57.5, -18.2715, 57.85, -18.2715};
              test_clip2D(2,3,4,p,4,q);
              std::vector<polygon> differences;
              weiler_clip_difference(2,2,4,p,4,q, differences);
            } 
            { // One edge of the cut polygon overlaps with the subject polygon internally
              double p[8] = {58, -19.8, 56, -19.8, 56, -18, 58, -18 };
              double q[8] = {58, -18.65, 57.5, -18.65, 57.5, -18.2715, 58, -18.2715};
              test_clip2D(2,3,4,p,4,q);
            } 
            {
              double p[6] = {-55.5,  -38.0685, -55,  -38.3136, -55,  -36.7609};
              double q[8] = { -55.5,  -38.5, -55,  -38.5, -55,  -38, -55.5,  -38};
              test_clip2D(2,3,3,p,4,q);
            }
            { // contains conincidental line segment
              double p[6] = {50.2073, 3.60505, 49,  3, 49.5,  3};
              double q[8] = {49,  3, 49.5,  3, 49.5,  3.25, 49,  3.25};
              test_clip2D(2,3,3,p,4,q);
            } 
            { // s contains c with zero common vertex
              double p[6] = {38.3779,  55.0379, 41.0815,  54.3093, 39.3743,  52.9098};
              double q[8] = {40,  53.5, 39.5,  53.5, 39.5,  54, 40,  54};
              test_clip2D(2,3,3,p,4,q); 
            }
  }
  return 0;
}
