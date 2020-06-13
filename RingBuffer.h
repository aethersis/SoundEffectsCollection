/*
  ==============================================================================

    RingBuffer.h
    Created: 23 May 2020 12:32:09pm
    Author:  mw

  ==============================================================================
*/

#pragma once
#include <cstdint>
#include <memory>
#include <cmath>

class RingBuffer
{
public:
    RingBuffer(size_t max_length) :
        max_length(max_length),
        current_length(max_length),
        buffer(new float[max_length]),
        head(max_length - 1u),
        tail(0)
    {
        for (size_t i = 0; i < max_length; ++i) {
            buffer[i] = 0.f;
        }
    }

    bool resize(size_t new_length)
    {
        if (new_length > max_length) {
            return false;
        }

        current_length = new_length;
        return true;
    }

    void update(float sample)
    {
        buffer[head] = sample;
        head = (head + 1u) % current_length;
        tail = (tail + 1u) % current_length;
    }

    float get(float offset = 0u)
    {
        // Yep, fractional length is possible :) This is useful for variations of Karplus Strong model (keytar) algorithms, time-stretching based pitch shift etc.
        float integerPart;
        float fracPart = modf(offset, &integerPart);

        return (1.f - fracPart) * buffer[(tail + (uint8_t)offset) % current_length] + 
               (      fracPart) * buffer[(tail + (uint8_t)offset + 1u) % current_length];
    }

private:
    float* buffer;
    size_t max_length;
    size_t current_length;
    size_t head;
    size_t tail;
};