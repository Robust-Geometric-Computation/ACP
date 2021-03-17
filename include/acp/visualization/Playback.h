#ifndef PLAYBACK_H
#define PLAYBACK_H

#include "acp/visualization/Playback.h"

#include <cstdio>
#include <ctime>
#include <iomanip>
#include <map>
#include <set>
#include <unordered_set>
#include <vector>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <utility>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <unistd.h>
#include <glm/gtc/quaternion.hpp> 
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/matrix_decompose.hpp>


const int ARRAY_SIZE = 16;


class Playback{
    

    

    public:
        Playback(const char*);
        bool isEmpty = true;
        bool started = false;
        bool interpolating = false;
        bool keep_interpolating = false;
        int pose_index;
        int pose_amount;
        float counter;
        float steps;
        
        bool isRegrid;
        bool update_model = false;

        glm::mat4 M1;
        glm::mat4 M2;
        glm::mat4 M_inverse;
        glm::mat4 identity_matrix = glm::mat4(1.f);

        glm::mat4 getFirstPose();
        void setTransformation();

        void updateM1(glm::mat4);
        void updateM2(glm::mat4);

        glm::mat4 getM2();
        glm::mat4 getS();
        glm::mat4 getR();
        glm::mat4 getT();

        float getCounter();
        float getSteps();
        

        glm::mat4 extracted_S;
        glm::mat4 extracted_R;
        glm::mat4 extracted_T;
        void printMat(glm::mat4);
        void recordCurrentPose(const char*, glm::mat4);
        void interpolate(glm::mat4&, glm::mat4&, glm::mat4&);

        float getScaleFactor(glm::mat4);
        glm::mat4 matpow(glm::mat4, float);
        void make_trs();
        glm::mat4 trs = glm::mat4(1.f);
        glm::mat4 trsinv = glm::mat4(1.f);
        glm::mat4 S;
        glm::mat4 T;
        glm::mat4 R;

        glm::mat4 getTRS();
        glm::mat4 getTRSinv();
        glm::mat4 getM1();
    
        std::vector<bool> regrid;
        bool isIdentity(glm::mat4 pose);

    private:
        std::vector<glm::mat4> poses;
        std::vector<float> numsteps;  
        void readPose(const char*);
        glm::mat4 makeCpy(glm::mat4);
        
    
};

#endif 
