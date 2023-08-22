#pragma once

#include "torrey.h"
#include "vector.h"
#include "image.h"
#include <variant>
#include "flexception.h"
#include "parse_scene.h"

struct ConstantTexture {
    Vector3 color;
};

struct ImageTexture {
    Image3 image;
    Real uscale = 1, vscale = 1;
    Real uoffset = 0, voffset = 0;
};

using Texture = std::variant<ConstantTexture, ImageTexture>;

inline Vector3 eval(const Texture &texture, const Vector2 &uv) {
    if (auto *constant = std::get_if<ConstantTexture>(&texture)) {
        return constant->color;
    } else if (auto *image = std::get_if<ImageTexture>(&texture)) {
        const Image3 &img = image->image;
        Real u = modulo(image->uscale * uv[0] + image->uoffset, Real(1)) *
            img.width;
        Real v = modulo(image->vscale * uv[1] + image->voffset, Real(1)) *
            img.height;
        int ufi = modulo(int(u), img.width);
        int vfi = modulo(int(v), img.height);
        int uci = modulo(ufi + 1, img.width);
        int vci = modulo(vfi + 1, img.height);
        Real u_off = u - ufi;
        Real v_off = v - vfi;
        Vector3 val_ff = img(ufi, vfi);
        Vector3 val_fc = img(ufi, vci);
        Vector3 val_cf = img(uci, vfi);
        Vector3 val_cc = img(uci, vci);
        return val_ff * (1 - u_off) * (1 - v_off) +
               val_fc * (1 - u_off) *      v_off +
               val_cf *      u_off  * (1 - v_off) +
               val_cc *      u_off  *      v_off;
    } else {
        Error("Unhandled Texture type");
        return Vector3{0, 0, 0};
    }
}

inline Texture parsed_color_to_texture(const ParsedColor &color) {
    if (auto *constant = std::get_if<Vector3>(&color)) {
        return ConstantTexture{*constant};
    } else if (auto *image = std::get_if<ParsedImageTexture>(&color)) {
        Image3 img = imread3(image->filename);
        return ImageTexture{img,
            image->uscale, image->vscale,
            image->uoffset, image->voffset};
    } else {
        Error("Unhandled ParsedColor");
    }
}