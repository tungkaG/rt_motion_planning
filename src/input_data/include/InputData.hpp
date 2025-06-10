#ifndef INPUTDATA_HPP
#define INPUTDATA_HPP

#include <iostream>


class InputData
{
public:

    InputData();

    InputData(const float& s_,
                    const float& ss_,
                    const float& sss_,
                    const float& d_,
                    const float& dd_,
                    const float& ddd_,
                    const float& velocity_,
                    const float& timestep_,
                    const float& desired_velocity_);

    float s;                   
    float ss;                   
    float sss;               
    float d;            
    float dd;        
    float ddd;

    float velocity;
    float timestep;
    float desired_velocity;
    
    void print(std::ostream& os) const;
};


#endif //InputData
