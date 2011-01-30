#ifndef dSFMT_wrapper_H
#define dSFMT_wrapper_H

// See notes in dSFMT_wrapper.cpp.

extern "C"
{

    void init_gen_rand_(uint32_t &seed);

    double genrand_close_open_(void);

    void fill_array_close_open_(double array[], int &size);

}

#endif // dSFMT_wrapper_H
