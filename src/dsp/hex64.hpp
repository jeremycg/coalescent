#pragma once

#include <cstddef>
#include <cstdint>

namespace coalescent {

// Lossless fixed-width JSON representation for unsigned 64-bit state. Parsing
// is transactional and accepts either case; formatting is canonical lowercase.
class Hex64Codec {
public:
    enum {
        DIGITS = 16,
        TEXT_SIZE = DIGITS + 1
    };

    static void format(std::uint64_t value, char (&text)[TEXT_SIZE]) {
        static const char digits[] = "0123456789abcdef";
        for (int index = DIGITS - 1; index >= 0; --index) {
            text[index] = digits[value & UINT64_C(0xf)];
            value >>= 4u;
        }
        text[DIGITS] = '\0';
    }

    static bool parse(const char* text, std::size_t length,
                      std::uint64_t& destination) {
        if (!text || length != static_cast<std::size_t>(DIGITS))
            return false;

        std::uint64_t parsed = 0u;
        for (int index = 0; index < DIGITS; ++index) {
            const char character = text[index];
            unsigned digit = 0u;
            if (character >= '0' && character <= '9')
                digit = static_cast<unsigned>(character - '0');
            else if (character >= 'a' && character <= 'f')
                digit = static_cast<unsigned>(character - 'a' + 10);
            else if (character >= 'A' && character <= 'F')
                digit = static_cast<unsigned>(character - 'A' + 10);
            else
                return false;
            parsed = (parsed << 4u) | static_cast<std::uint64_t>(digit);
        }

        destination = parsed;
        return true;
    }

    static bool parseCString(const char* text, std::uint64_t& destination) {
        if (!text)
            return false;
        std::size_t length = 0u;
        while (length <= static_cast<std::size_t>(DIGITS) && text[length] != '\0')
            ++length;
        return parse(text, length, destination);
    }
};

} // namespace coalescent
