//
// Created by Allison Mason on 2/23/24 from GNPS_Scoring/SquareRootNormalize.
//

#ifndef SCALENORMALIZE_H
#define SCALENORMALIZE_H
#include <vector>


class ScaleNormalize {
public:
    inline std::vector<double> Normalize(const std::vector<double> value)
    {
        
        const size_t count = value.size();
        std::vector<double> new_values(value.size());
        double max = 0;
 
        for (size_t i = 0; i < count; i++) {
            if (value[i] > max) {
                max = value[i];
            }
        }

        for(size_t i = 0; i < count; i++)
        {
            new_values[i] = (value[i] / max);
        }
        return new_values;
    }
};



#endif // SCALENORMALIZE_H
