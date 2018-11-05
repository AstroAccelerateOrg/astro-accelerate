//
//  aa_ddtr_plan.hpp
//  aapipeline
//
//  Created by Cees Carels on Tuesday 23/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef aa_ddtr_plan_hpp
#define aa_ddtr_plan_hpp

#include <stdio.h>
#include <vector>

class aa_ddtr_plan {
public:
    aa_ddtr_plan() {
        
    }
    
    struct dm {
        float low;
        float high;
        float step;
        int inBin;
        int outBin;
    };
    
    const bool add_dm(const float &low, const float &high, const float &step, const int &inBin, const int &outBin) {
        const dm tmp = {low, high, step, inBin, outBin};
        m_user_dm.push_back(std::move(tmp));
        return true;
    }
    
    const size_t range() const {
        return m_user_dm.size();
    }
    
    const dm user_dm(const size_t &i) const {
        return m_user_dm.at(i);
    }
    
private:
    std::vector<dm> m_user_dm;
    
};

#endif /* aa_ddtr_plan_hpp */
