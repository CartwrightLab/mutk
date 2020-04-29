/*
# Copyright (c) 2020 Reed A. Cartwright <reed@cartwright.ht>
#
# This file is part of the Ultimate Source Code Project.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
*/

#ifndef MUTK_UTILITY_HPP
#define MUTK_UTILITY_HPP

#include <string>
#include <filesystem>
#include <iostream>
#include <fstream>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/tokenizer.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/distance.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/as_literal.hpp>
#include <boost/range/reference.hpp>
#include <boost/range/size.hpp>

namespace mutk {
namespace utility {


// Key search functions
template<typename Range, typename Value>
inline
size_t find_position(const Range& r, const Value& v) {
    return boost::distance(boost::find<boost::return_begin_found>(r, v));
}

template<typename Range, typename Pred>
inline
size_t find_position_if(const Range& r, Pred pred) {
    return boost::distance(boost::find_if<boost::return_begin_found>(r, pred));
}

template<typename A, typename B>
std::size_t key_switch(const A &ss, const B &keys) {
    std::size_t ret = find_position_if(keys,
        [&](typename boost::range_reference<const B>::type k) {
            return boost::algorithm::istarts_with(k, ss);
        });
    return (ret != boost::size(keys)) ? ret : static_cast<std::size_t>(-1);
}

template<typename A, typename B>
std::size_t key_switch_iequals(const A &ss, const B &keys) {
    std::size_t ret = find_position_if(keys,
        [&](typename boost::range_reference<const B>::type k) {
            return boost::algorithm::iequals(k, ss);
        });
    return (ret != boost::size(keys)) ? ret : static_cast<std::size_t>(-1);
}

template<typename A, typename B>
typename boost::range_reference<const B>::type
key_switch_tuple(const A &ss, const B &keys,
    typename boost::range_reference<const B>::type &default_value) {
    auto it = boost::find_if<boost::return_found>(keys,
        [&](typename boost::range_reference<const B>::type k) {
            return boost::algorithm::istarts_with(std::get<0>(k), ss);
        });
    return (it != boost::end(keys)) ? *it : default_value;
}


// Tokenization Functions
namespace detail {
    using token_function = boost::char_separator<char>;
    template<typename Range>
    using char_tokenizer = boost::tokenizer<token_function, typename boost::range_iterator<const Range>::type>;
};

template<typename Range>
inline
detail::char_tokenizer<Range>
make_tokenizer(const Range& text, const char *sep = "\t", const char *eol = "\n") {
    detail::token_function f(sep, eol, boost::keep_empty_tokens);
    return {boost::as_literal(text),f};
}

template<typename Range>
inline
detail::char_tokenizer<Range>
make_tokenizer_dropempty(const Range& text, const char *sep = "\t", const char *eol = "\n") {
    detail::token_function f(sep, eol, boost::drop_empty_tokens);
    return {boost::as_literal(text),f};
}

// decode a percent encoded string
namespace detail {
void percent_decode_core(std::string *str, size_t start);
} // namespace detail

inline
std::string percent_decode(std::string str) {
    auto pos = str.find('%');
    if(pos != std::string::npos) {
        detail::percent_decode_core(&str, pos);
    }
    return std::move(str);
}

// extracts extension and filename from both file.foo and ext:file.foo
struct file_type_t {
    std::string path;
    std::string type_ext;
    std::string compress_ext;
};

// returns {ext, file.foo}
// trims whitespace as well
file_type_t extract_file_type(std::string path);

// wrapper for input files that might be either a file or stdin
class File {
public:
    File() = default;
    
    File(File&& other);

    File& operator=(File&& other);

    explicit File(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in) {
        Open(filename, mode);
    }

    bool Open(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in);

    bool Attach(std::streambuf *buffer) {
        stream_.rdbuf(buffer);
        stream_.unsetf(std::ios_base::skipws);
        return is_open();
    }

    bool is_open() const { return stream_.rdbuf() != nullptr; }
    operator bool() const { return is_open(); }

    std::string path() const { return path_.native(); }

protected:
    std::iostream stream_{nullptr};
    std::filesystem::path path_;
    std::string type_label_;
    std::string compression_;

private:
    std::unique_ptr<std::streambuf> buffer_;    
};

inline
File::File(File&& other) : 
    buffer_{std::move(other.buffer_)},
    path_{std::move(other.path_)},
    type_label_{std::move(other.type_label_)}
{
    std::streambuf *buffer = other.stream_.rdbuf();
    Attach(buffer);
    other.Attach(nullptr);
} 

inline
File& File::operator=(File&& other) {
    if(this == &other) {
        return *this;
    }

    buffer_ = std::move(other.buffer_);
    path_ = std::move(other.path_);
    type_label_ = std::move(other.type_label_);

    std::streambuf *buffer = other.stream_.rdbuf();
    Attach(buffer);
    other.Attach(nullptr);

    return *this;
}

class BinaryFile : public File {
public:
    BinaryFile() = default;

    BinaryFile(BinaryFile&&) = default;
    BinaryFile& operator=(BinaryFile&&) = default;
    
    explicit BinaryFile(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in) {
        Open(filename, mode);
    }

    void Open(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in) {
        File::Open(filename, mode | std::ios::binary);
    }
};

template<typename X, typename T = std::char_traits<X>, typename A = std::allocator<X>>
inline
std::optional<std::basic_string<X,T,A>> slurp(std::basic_ifstream<X,T>& input) {
    std::basic_string<X,T,A> ret;
    if(!input) {
        return std::nullopt;
    }
    input.seekg(0, std::ios::end);
    ret.resize(input.tellg());
    input.seekg(0, std::ios::beg);
    input.read(&ret[0], ret.size());
    return ret;
}

inline std::optional<std::string> slurp(const std::filesystem::path& path, std::ios_base::openmode mode = std::ios_base::in) {
    std::ifstream in{path,mode};
    return slurp(in);
}

} // namespace utility
} // namespace mutk

#endif // MUTK_UTILITY_HPP