#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "AstroAccelerate/device_peak_find.h"

#include "managed_allocator.hpp"
#include <vector>
#include <algorithm>
#include <iostream>

template <typename T>
void debug_print_output(T & input, int width, int height) {
    for(int i=0; i < height; ++i) {
        if (i > 0 && i % 32 == 0) { std::cout << "  "; for (int j=0; j < width + ((width - 1) / 32); ++j) { std::cout << " -"; }  std::cout << "\n";} 
        std::cout.width(2);
        std::cout << i; 
        std::cout.width(0);
        for(int j=0; j < width; ++j) {
            if (j > 0 && j % 32 == 0) std::cout << " |";
            std::cout << " " << input[i*width+j];
        }
        std::cout << "\n";
    }
    std::cout << std::flush;
}

template<typename Func>
void test_peak_at_loc(std::vector<float, managed_allocator<float>> & input, std::vector<unsigned short, managed_allocator<unsigned short>> & output, int peak_loc, Func peakfind_func)
{
  input.resize(9);
  output.resize(9);
  for(int i=0; i < 9; ++i) {
    input[i] = (i==peak_loc) ? 100.0f : 0.0f;
  }
  peakfind_func(input.data(), 0, 3, 3, output.data());
  for (int i=0; i < 9; ++i) {
    auto val = (i==peak_loc) ? 1u : 0u;
    REQUIRE(output[i] == val);
  }
}

TEST_CASE("Peak Find v1", "[peakfind]") {
  std::vector<float, managed_allocator<float>> data;
  std::vector<unsigned short, managed_allocator<unsigned short>> output;

  SECTION ("Empty data is ok") {
    peakfind(data.data(), 0, 0, 0, output.data());
  }

  SECTION ("Single element is peak") {
    data.resize(1);
    output.resize(1);
    data[0] = 1.0;
    peakfind(data.data(), 0, 1, 1, output.data());

    REQUIRE(output[0] == 1);
  }

  SECTION ("Peak in top left") {
    test_peak_at_loc(data, output, 0, peakfind);
  }

  SECTION ("Peak on top edge") {
    test_peak_at_loc(data, output, 1, peakfind);
  }
  SECTION ("Peak in top right") {
    test_peak_at_loc(data, output, 2, peakfind);
  }
  SECTION ("Peak on left edge") {
    test_peak_at_loc(data, output, 3, peakfind);
  }
  SECTION ("Peak in centre") {
    test_peak_at_loc(data, output, 4, peakfind);
  }
  SECTION ("Peak on right edge") {
    test_peak_at_loc(data, output, 5, peakfind);
  }
  SECTION ("Peak in bottom left") {
    test_peak_at_loc(data, output, 6, peakfind);
  }
  SECTION ("Peak on bottom edge") {
    test_peak_at_loc(data, output, 7, peakfind);
  }
  SECTION ("Peak in bottom right") {
    test_peak_at_loc(data, output, 8, peakfind);
  }

  SECTION ("Peak in neighbouring block") {
    data.resize(1024*2*3);
    output.resize(data.size());
    std::fill(data.begin(), data.end(), 0.0f);

    SECTION("Right BC") {
        data[1024] = 10000.0f;
        data[1025] = 1.0f;
        peakfind(data.data(), 0, 1024*2, 3, output.data());
        for(std::size_t i=0; i < data.size(); ++i) {
            auto val = (i==1024) ? 1u: 0u;
            REQUIRE(output[i]==val);
        }
    } 
    SECTION("Left BC") {
        data[1024] = 1.0f;
        data[1024+1] = 10000.0f;
        peakfind(data.data(), 0, 1024*2, 3, output.data());
        for(std::size_t i=0; i < data.size(); ++i) {
            auto val = (i==1025) ? 1u: 0u;
            REQUIRE(output[i]==val);
        }
    }
  }
  SECTION ("Peak across block corners") {
     data.resize(1024 *4);
     output.resize(data.size());
     std::fill(data.begin(), data.end(), 0.0f);

     SECTION("Top left block BR corner") {
         data[64*32-32-1] = 100.0f;
         data[64*32+32] = 1.0f;
         peakfind(data.data(), 0, 64, 64, output.data());
         for(std::size_t i=0; i < data.size(); ++i) {
            auto val = (i==64*32-32-1) ? 1u: 0u;
            REQUIRE(output[i]==val);
         }
     }
     SECTION("Top right block BL corner") {
         data[64*32-32] = 100.0f;
         data[64*32+32-1] = 1.0f;
         peakfind(data.data(), 0, 64, 64, output.data());
         for(std::size_t i=0; i < data.size(); ++i) {
            auto val = (i==64*32-32) ? 1u: 0u;
            REQUIRE(output[i]==val);
         }
     }
     SECTION("Bottom right block TL corner") {
         data[64*32-32-1] = 1.0f;
         data[64*32+32] = 100.0f;
         peakfind(data.data(), 0, 64, 64, output.data());
         for(std::size_t i=0; i < data.size(); ++i) {
            auto val = (i==64*32+32) ? 1u: 0u;
            REQUIRE(output[i]==val);
         }
     }
     SECTION("Bottom left block TR corner") {
         data[64*32-32] = 1.0f;
         data[64*32+32-1] = 100.0f;
         peakfind(data.data(), 0, 64, 64, output.data());
         for(std::size_t i=0; i < data.size(); ++i) {
            auto val = (i==64*32+32-1) ? 1u: 0u;
            REQUIRE(output[i]==val);
         }
     }
    }
}

TEST_CASE("Peak Find v2", "[peakfind]") {
  std::vector<float, managed_allocator<float>> data;
  std::vector<unsigned short, managed_allocator<unsigned short>> output;

  SECTION ("Empty data is ok") {
    peakfind_v2(data.data(), 0, 0, 0, output.data());
  }

  SECTION ("Single element is peak") {
    data.resize(1);
    output.resize(1);
    data[0] = 1.0;
    peakfind_v2(data.data(), 0, 1, 1, output.data());

    REQUIRE(output[0] == 1);
  }

  SECTION ("Peak in top left") {
    test_peak_at_loc(data, output, 0, peakfind_v2);
  }

  SECTION ("Peak on top edge") {
    test_peak_at_loc(data, output, 1, peakfind_v2);
  }
  SECTION ("Peak in top right") {
    test_peak_at_loc(data, output, 2, peakfind_v2);
  }
  SECTION ("Peak on left edge") {
    test_peak_at_loc(data, output, 3, peakfind_v2);
  }
  SECTION ("Peak in centre") {
    test_peak_at_loc(data, output, 4, peakfind_v2);
  }
  SECTION ("Peak on right edge") {
    test_peak_at_loc(data, output, 5, peakfind_v2);
  }
  SECTION ("Peak in bottom left") {
    test_peak_at_loc(data, output, 6, peakfind_v2);
  }
  SECTION ("Peak on bottom edge") {
    test_peak_at_loc(data, output, 7, peakfind_v2);
  }
  SECTION ("Peak in bottom right") {
    test_peak_at_loc(data, output, 8, peakfind_v2);
  }

  SECTION ("Peak in neighbouring block") {
    data.resize(1024*2*3);
    output.resize(data.size());
    std::fill(data.begin(), data.end(), 0.0f);

    SECTION("Right BC") {
        data[1024] = 10000.0f;
        data[1025] = 1.0f;
        peakfind_v2(data.data(), 0, 1024*2, 3, output.data());
        for(std::size_t i=0; i < data.size(); ++i) {
            auto val = (i==1024) ? 1u: 0u;
            REQUIRE(output[i]==val);
        }
    } 
    SECTION("Left BC") {
        data[1024] = 1.0f;
        data[1024+1] = 10000.0f;
        peakfind_v2(data.data(), 0, 1024*2, 3, output.data());
        for(std::size_t i=0; i < data.size(); ++i) {
            auto val = (i==1025) ? 1u: 0u;
            REQUIRE(output[i]==val);
        }
    }
  }
  SECTION ("Peak across block corners") {
     data.resize(1024 *4);
     output.resize(data.size());
     std::fill(data.begin(), data.end(), 0.0f);

     SECTION("Top left block BR corner") {
         data[64*32-32-1] = 100.0f;
         data[64*32+32] = 1.0f;
         peakfind_v2(data.data(), 0, 64, 64, output.data());
         for(std::size_t i=0; i < data.size(); ++i) {
            auto val = (i==64*32-32-1) ? 1u: 0u;
            REQUIRE(output[i]==val);
         }
     }
     SECTION("Top right block BL corner") {
         data[64*32-32] = 100.0f;
         data[64*32+32-1] = 1.0f;
         peakfind_v2(data.data(), 0, 64, 64, output.data());
         for(std::size_t i=0; i < data.size(); ++i) {
            auto val = (i==64*32-32) ? 1u: 0u;
            REQUIRE(output[i]==val);
         }
     }
     SECTION("Bottom right block TL corner") {
         data[64*32-32-1] = 1.0f;
         data[64*32+32] = 100.0f;
         peakfind_v2(data.data(), 0, 64, 64, output.data());
         for(std::size_t i=0; i < data.size(); ++i) {
            auto val = (i==64*32+32) ? 1u: 0u;
            REQUIRE(output[i]==val);
         }
     }
     SECTION("Bottom left block TR corner") {
         data[64*32-32] = 1.0f;
         data[64*32+32-1] = 100.0f;
         peakfind_v2(data.data(), 0, 64, 64, output.data());
         for(std::size_t i=0; i < data.size(); ++i) {
            auto val = (i==64*32+32-1) ? 1u: 0u;
            REQUIRE(output[i]==val);
         }
     }
    }
}

TEST_CASE("Peak Find v3", "[peakfind]") {
  std::vector<float, managed_allocator<float>> data;
  std::vector<unsigned short, managed_allocator<unsigned short>> output;

  SECTION ("Empty data is ok") {
    peakfind_v3(data.data(), 0, 0, 0, output.data());
  }

  SECTION ("Single element is peak") {
    data.resize(1);
    output.resize(1);
    data[0] = 1.0;
    peakfind_v3(data.data(), 0, 1, 1, output.data());

    REQUIRE(output[0] == 1);
  }

  SECTION ("Peak in top left") {
    test_peak_at_loc(data, output, 0, peakfind_v3);
  }

  SECTION ("Peak on top edge") {
    test_peak_at_loc(data, output, 1, peakfind_v3);
  }
  SECTION ("Peak in top right") {
    test_peak_at_loc(data, output, 2, peakfind_v3);
  }
  SECTION ("Peak on left edge") {
    test_peak_at_loc(data, output, 3, peakfind_v3);
  }
  SECTION ("Peak in centre") {
    test_peak_at_loc(data, output, 4, peakfind_v3);
  }
  SECTION ("Peak on right edge") {
    test_peak_at_loc(data, output, 5, peakfind_v3);
  }
  SECTION ("Peak in bottom left") {
    test_peak_at_loc(data, output, 6, peakfind_v3);
  }
  SECTION ("Peak on bottom edge") {
    test_peak_at_loc(data, output, 7, peakfind_v3);
  }
  SECTION ("Peak in bottom right") {
    test_peak_at_loc(data, output, 8, peakfind_v3);
  }

  SECTION ("Peak in neighbouring block") {
    data.resize(1024*2*3);
    output.resize(data.size());
    std::fill(data.begin(), data.end(), 0.0f);

    SECTION("Right BC") {
        data[1024] = 10000.0f;
        data[1025] = 1.0f;
        peakfind_v3(data.data(), 0, 1024*2, 3, output.data());
        for(std::size_t i=0; i < data.size(); ++i) {
            auto val = (i==1024) ? 1u: 0u;
            REQUIRE(output[i]==val);
        }
    } 
    SECTION("Left BC") {
        data[1024] = 1.0f;
        data[1024+1] = 10000.0f;
        peakfind_v3(data.data(), 0, 1024*2, 3, output.data());
        for(std::size_t i=0; i < data.size(); ++i) {
            auto val = (i==1025) ? 1u: 0u;
            REQUIRE(output[i]==val);
        }
    }
  }
  SECTION ("Peak across block corners") {
     data.resize(1024 *4);
     output.resize(data.size());
     std::fill(data.begin(), data.end(), 0.0f);

     SECTION("Top left block BR corner") {
         data[64*32-32-1] = 100.0f;
         data[64*32+32] = 1.0f;
         peakfind_v3(data.data(), 0, 64, 64, output.data());
         for(std::size_t i=0; i < data.size(); ++i) {
            auto val = (i==64*32-32-1) ? 1u: 0u;
            REQUIRE(output[i]==val);
         }
     }
     SECTION("Top right block BL corner") {
         data[64*32-32] = 100.0f;
         data[64*32+32-1] = 1.0f;
         peakfind_v3(data.data(), 0, 64, 64, output.data());
         for(std::size_t i=0; i < data.size(); ++i) {
            auto val = (i==64*32-32) ? 1u: 0u;
            REQUIRE(output[i]==val);
         }
     }
     SECTION("Bottom right block TL corner") {
         data[64*32-32-1] = 1.0f;
         data[64*32+32] = 100.0f;
         peakfind_v3(data.data(), 0, 64, 64, output.data());
         for(std::size_t i=0; i < data.size(); ++i) {
            auto val = (i==64*32+32) ? 1u: 0u;
            REQUIRE(output[i]==val);
         }
     }
     SECTION("Bottom left block TR corner") {
         data[64*32-32] = 1.0f;
         data[64*32+32-1] = 100.0f;
         peakfind_v3(data.data(), 0, 64, 64, output.data());
         for(std::size_t i=0; i < data.size(); ++i) {
            auto val = (i==64*32+32-1) ? 1u: 0u;
            REQUIRE(output[i]==val);
         }
     }
    }
}
