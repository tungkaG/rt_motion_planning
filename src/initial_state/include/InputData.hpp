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
                    const float& ddd_);

    float s;                   
    float ss;                   
    float sss;               
    float d;            
    float dd;        
    float ddd;

    x_0.velocity and x_0.timestep, and desired_velocity (see what are required to calculate d_sampling, t_sampling and v_sampling)
    
    void print(std::ostream& os) const;
};


#endif //InputData
