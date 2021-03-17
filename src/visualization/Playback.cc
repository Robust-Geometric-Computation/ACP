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

#include <tgmath.h> 

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


const int ARRAY_SIZE = 16;


class Playback{
    

    public:
        Playback(const char*);
        bool isEmpty = true;
        bool started = false;
        bool interpolating = false;
        int pose_index;
        int pose_amount;
        float counter;
        float steps;
        bool isRegrid;
        bool update_model = false;
        bool keep_interpolating = false;

        glm::mat4 M1;
        glm::mat4 M2;
        glm::mat4 M_inverse;
        glm::mat4 identity_matrix = glm::mat4(1.f);

        glm::mat4 getFirstPose();
        void setTransformation();

        void updateM1(glm::mat4);
        void updateM2(glm::mat4);

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
        bool isIdentity(glm::mat4 pose);

        std::vector<bool> regrid;

    private:
        std::vector<glm::mat4> poses;
        std::vector<float> numsteps;  
        void readPose(const char*);
        glm::mat4 makeCpy(glm::mat4);
        
        
    
};

Playback::Playback(const char* filename)
{
    readPose(filename);
    this->isEmpty = false;
    this->started = false;
    pose_index = 1;
}

glm::mat4 Playback::getFirstPose()
{   
    this->started = true;
    return makeCpy(poses[0]);
}


void Playback::readPose(const char* filename)
{
    std::ifstream pose_reader;
    

	pose_reader.open(filename);
	
    if(pose_reader.fail())
    {
        std::cout << "Could not open file: \"" << filename << "\". Exiting.\n\n";
        exit(0);
    }
	
    if (!pose_reader)
    {
      std::cout << "Error: Can't open the file. \n";
      exit(1);
    }

    float num;

    float file_matrix[ARRAY_SIZE];

    int indexing = 0;

    while(pose_reader >> num)
    {
        if(indexing == ARRAY_SIZE)
        {
            glm:: mat4 pose_matrix = glm::make_mat4(file_matrix);
            
            pose_matrix = glm::transpose(pose_matrix);

            poses.push_back(pose_matrix);
            if(num < 0)
            {
                regrid.push_back(true);
                pose_reader >> num;
                numsteps.push_back(num);
            }
            else
            {
                regrid.push_back(false);
                numsteps.push_back(num);
            }
            indexing = -1;
        }
        else
        {
            file_matrix[indexing] = num;
        }
        
        indexing++;
    } 
    pose_reader.close();
    pose_amount = poses.size();
}

void Playback::setTransformation()
{   
    steps = numsteps[pose_index]; 

    isRegrid = regrid[pose_index];
    
    updateM2(poses[pose_index]);

    glm::mat4 M3;
    glm::fquat R_Quat;
    glm::fquat identity_quat = glm::quat_cast(glm::mat4(1.f));

    float scale_factor;
    float nth_root;
    M3 = M2 * glm::inverse(M1);

    scale_factor = getScaleFactor(M3);
   
    if(std::abs(scale_factor - 1.0f) < 0.0001)
    {
        scale_factor = 1.0f;
    }

   
    nth_root = pow(scale_factor, 1.0f/steps);
    S = glm::scale(glm::mat4(1.f), glm::vec3(scale_factor, scale_factor, scale_factor));
    extracted_S = glm::scale(glm::mat4(1.f), glm::vec3(nth_root, nth_root, nth_root));

    R = {
            {M3[0][0] / scale_factor, M3[0][1] / scale_factor, M3[0][2] /scale_factor, 0.0f},
            {M3[1][0] / scale_factor, M3[1][1] / scale_factor, M3[1][2] / scale_factor, 0.0f},
            {M3[2][0] / scale_factor, M3[2][1] / scale_factor, M3[2][2] / scale_factor, 0.0f},
            {0.0f, 0.0f, 0.0f, 1.0f}
        };


    R_Quat = glm::quat_cast(R);

    R_Quat = glm::slerp(identity_quat, R_Quat, 1.0f / steps);

    R_Quat = glm::normalize(R_Quat);

    extracted_R = glm::mat4_cast(R_Quat);

    T = glm::translate(glm::mat4(1.f), glm::vec3(M3[3][0], M3[3][1], M3[3][2]));

    glm::vec3 t = glm::vec3(M3[3][0] / steps, M3[3][1] / steps, M3[3][2] / steps);

    extracted_T = glm::translate(glm::mat4(1.f), t);

    make_trs();

    counter = 0;

    pose_index++;

    if(pose_index == pose_amount)
    {
        std::cout << "\nReached end of file.\n\n";
        
        this->isEmpty = true;

       this->keep_interpolating = false;
        return;
    }
}

float Playback::getScaleFactor(glm::mat4 pose)
{
  glm::vec3 col1 = {pose[0][0], pose[1][0], pose[2][0]};
  glm::vec3 col2 = {pose[0][1], pose[1][1], pose[2][1]};
  glm::vec3 col3 = {pose[0][2], pose[1][2], pose[2][2]};

  float average = (glm::length(col1) + glm::length(col2) + glm::length(col3)) / 3.0;

  return average;
}



void Playback::updateM1(glm::mat4 m1)
{
    M1 = m1;
    return;
}

void Playback::updateM2(glm::mat4 m2)
{
    M2 = m2;
    return;
}
glm::mat4 Playback::getM1()
{
    return makeCpy(M1);
}


glm::mat4 Playback::makeCpy(glm::mat4 copmat)
{
    float arr[16];
    int x = 0;
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            arr[x] = copmat[i][j];
            x++;
        }
    }

    return glm::make_mat4(arr);
}

void Playback::printMat(glm::mat4 po)
{
         for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            std::cout << po[i][j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n\n";
}

glm::mat4 Playback::matpow(glm::mat4 pose, float power)
{   
    glm::mat4 idmat = glm::mat4(1.f);

    for(int i = 0; i < power; i++)
    {
        idmat = pose * idmat;
    }

    return idmat;
}

void Playback::make_trs()
{
  
    glm::mat4 t = glm::mat4(1.f);
    glm::mat4 s = extracted_S;
    glm::mat4 r = extracted_R;

    glm::mat4 RS = R * S;
    glm::mat4 rs = r * s;
    glm::mat4 RS_inv;

    //IF T IS IDENTITY GET OUT OF HERE

    float difx = std::abs(T[3][0]);
    float dify = std::abs(T[3][1]);
    float difz = std::abs(T[3][2]);

    if(difx < 0.001 && dify < 0.001 && difz < 0.001)
    {
        std::cout << "T is identity. returning \n\n";
        trs = r * s;
        trsinv = glm::transpose(glm::inverse(trs));
        return;
    }


    S = glm::mat4(1.f);

    for(int i = 0; i < this->steps; i++)
    {
        S = (rs * S) + glm::mat4(1.f);
    }

    S[3][3] = 1.0f;

    if(isIdentity(S) && isIdentity(R))
    {   
        t = extracted_T;
        trs = t;
        trsinv = glm::transpose(glm::inverse(trs));
    }
    else
    {   
        //t = RS_inv * rs * T * t;
        t = glm::inverse(S) * T;

        t = glm::translate(glm::mat4(1.f), glm::vec3(t[3][0], t[3][1], t[3][2]));
        trs = t * r * s;
        trsinv = glm::transpose(glm::inverse(trs));
    }
        

    
    //std::cout << "t : \n\n"; printMat(t);
    //std::cout << " trs : \n\n"; printMat(trs);
    //std::cout << "TRS * M4 : \n\n"; printMat(T * R * S * M1);
    //std::cout << "trs^n * M4 : \n\n"; printMat(matpow(trs, steps) * M1);
}

glm::mat4 Playback::getTRS()
{
    return makeCpy(trs);
}

glm::mat4 Playback::getTRSinv()
{
    return makeCpy(trsinv);
}

bool Playback::isIdentity(glm::mat4 pose)
{
    bool check = true;
    glm::mat4 identity = glm::mat4(1.f);

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            if(std::abs(pose[i][j] - identity[i][j]) > 0.00001)
            {
                check = false;
            }

        }
    }

    return check;
}

float Playback::getCounter()
{
    return counter;
}

float Playback::getSteps()
{
    return steps;
}

void Playback::recordCurrentPose(const char* filename, glm::mat4 model)
{

  const int DEFAULT_NUMSTEPS = 100.0f;
  std::ofstream record_pose;
  
  
  record_pose.open(filename, std::ios_base::app);

  std::cout << "We are recording a pose.\n\n";

  for(int k = 0; k < 4; k++)
  {
      for(int l = 0; l < 4; l++)
      {
        record_pose << model[l][k] << " ";
      }
      record_pose << "\n";
  }

  record_pose << "\n" << DEFAULT_NUMSTEPS << "\n\n\n";

  record_pose.close();
}


void Playback::interpolate(glm::mat4 &model, glm::mat4 &MLines, glm::mat4 &model_inverse)
{   

    if(started == false){

        model = getFirstPose();
        MLines = model;

        update_model = true;

        return;
    }
    else{
        
        if(interpolating == false){  
            
            M1 = identity_matrix * model;
            M_inverse = identity_matrix * model_inverse;
            
            setTransformation();
            interpolating = true;

        }
        else{
        
            if(counter > steps - 1){
                model = M2;
                MLines = M2;
                update_model = true;
                interpolating = false;
            }
            else{
                
                M1 = trs * M1;
                M_inverse = trsinv * M_inverse;
                
                model = identity_matrix * M1;
                MLines = model;

                model_inverse = identity_matrix * M_inverse;
                counter++;
            }
            
        }
        
        
    }
    
    
    
} 
#endif
