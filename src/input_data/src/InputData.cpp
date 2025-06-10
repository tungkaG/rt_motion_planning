#include "InputData.hpp"

#include <iomanip>

InputData::InputData()
    : s ()
    , ss ()
    , sss () 
    , d ()
    , dd ()
    , ddd ()
{

}

InputData::InputData(const float& s_,
         const float& ss_,
         const float& sss_,
         const float& d_,
         const float& dd_,
         const float& ddd_,
         const float& velocity_,
         const float& timestep_,
         const float& desired_velocity_)
    : s (s_)
    , ss (ss_)
    , sss (sss_)
    , d (d_)
    , dd (dd_)
    , ddd (ddd_)
    , velocity (velocity_)
    , timestep (timestep_)
    , desired_velocity (desired_velocity_)
{

}

void InputData::print(std::ostream& os) const
{
    int width = 15;
    os << std::fixed << std::setprecision(5);

    os << "InputData:" << std::endl
       << std::setw(width) << "s"
       << std::setw(width) << "ss"
       << std::setw(width) << "sss"
       << std::setw(width) << "d"
       << std::setw(width) << "dd"
       << std::setw(width) << "ddd"
       << std::setw(width) << "velocity"
       << std::setw(width) << "timestep"
       << std::setw(width) << "desired_velocity" << std::endl;

    os << std::setw(width) << s
       << std::setw(width) << ss
       << std::setw(width) << sss
       << std::setw(width) << d
       << std::setw(width) << dd
       << std::setw(width) << ddd
       << std::setw(width) << velocity
       << std::setw(width) << timestep
       << std::setw(width) << desired_velocity << std::endl;
}