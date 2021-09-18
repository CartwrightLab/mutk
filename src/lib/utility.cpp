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

#include <mutk/utility.hpp>

#include <boost/filesystem/convenience.hpp>
#include <boost/algorithm/string/trim.hpp>

void mutk::utility::detail::percent_decode_core(std::string &str, size_t start) {
    assert(str[start] == '%');

    auto hex_decode = [](char x) -> int {
         if('0' <= x && x <= '9') {
            return x-'0';
        }
        if('A' <= x && x <= 'F') {
            return x-'A'+10;
        }
        if('a' <= x && x <= 'f') {
            return x-'a'+10;
        }
        return -1;
    };

    auto p = str.begin()+start;
    auto q = p;
    int a, b;
    do {
        if(++p == str.end()) {
            break;
        }
        a = hex_decode(*p);
        if(++p == str.end()) {
            break;
        }
        b = hex_decode(*p);
        if(a != -1 && b != -1) {
            *q++ = a*16+b;
        }
        for(++p; p!=str.end(); ++p) {
            if(*p == '%') {
                break;
            }
            *q++ = *p;
        }
    } while(p != str.end());
    
    str.erase(q,str.end());
}

// extracts extension and filename from both file.ext and ext:file.foo
// trims whitespace as well
mutk::utility::file_type_t mutk::utility::extract_file_type(std::string path) {
    using boost::algorithm::trim;
    using boost::algorithm::ends_with;
    const auto &npos = std::string::npos;

    auto is_compressed = [](const std::string &path) -> std::string {
        for(auto && s : {".gz", ".gzip", ".bgz"}) {
            if(ends_with(path, s)) {
                return {s+1};
            }
        }
        return {};
    };

    // Remove whitespace
    trim(path);

    // Format ext:path ???
    auto colon = path.find_first_of(':');
    if(colon != npos) {
        auto ext = path.substr(0, colon);
        path.erase(0,colon+1);
        // Format ext.gz:path ???
        auto gz = is_compressed('.'+ext);
        if(!gz.empty()) {
            // Format gz:path
            if(gz.size() == ext.size()) {
                ext.erase();
            } else {
                auto gz_pos = ext.size() - (gz.size()+1);
                ext.erase(gz_pos);
            }
        }
        return {path, ext, gz};
    }
    // Format path.gz ???
    auto gz = is_compressed(path);
    auto gz_pos = path.size() - (gz.empty() ? 0 : gz.size()+1);
    auto ext_pos = path.find_last_of('.', gz_pos-1);
    if(ext_pos == npos || ext_pos == 0) {
        return {path, {}, gz};
    }
    auto ext = path.substr(ext_pos+1, gz_pos-(ext_pos+1));
    return {path, ext, gz};
}

bool mutk::utility::File::Open(const std::string &filename, std::ios_base::openmode mode) {
    auto file_type = utility::extract_file_type(filename);
    path_ = file_type.path;
    type_label_ = file_type.type_ext;
    compression_ = file_type.compress_ext;

    if(path_.empty()) {
        // don't do anything, just setup type
    } else if(path_ != "-") {
        // if path is not "-" open a file
        buffer_.reset(new std::filebuf);
        std::filebuf *p = static_cast<std::filebuf*>(buffer_.get());
        p->open(path_.c_str(), mode);
        // if file is open, associate it with the stream
        if(p->is_open()) {
            return Attach(buffer_.get());
        }
    } else if((mode & std::ios_base::in) == (mode & std::ios_base::out)) {
        // can't do anything if both are set or none are set
    } else if(mode & std::ios_base::in) {
        return Attach(std::cin.rdbuf());
    } else {
        return Attach(std::cout.rdbuf());
    }
    return Attach(nullptr);
}