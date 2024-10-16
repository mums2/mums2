#ifndef UTILS
#define UTILS

#include <string>
#include <vector>
#include <cctype>
#include <iostream>
#include <Rcpp.h>
// #include <algorithm>

class Utils {
    public: 
        Utils() = default; 
        ~Utils() = default;
        
        std::string stringToLower(std::string s) {
            
            std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); } 
                  );
            return s;
        }

        std::vector<int> my_grep(std::vector<std::string> x, std::string pattern) {
            int n = x.size();

            std::vector<int> out;
            
            for (int i = 0; i < n; i++) {
 
                if (stringToLower(x[i]) != pattern) {
                    continue;
                }
                
                out.push_back(i);
            } 
            
            return out;
 
        }

};

#endif //UTILS