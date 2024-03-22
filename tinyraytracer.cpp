#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <time.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "model.h"
#include "geometry.h"

int envmap_width, envmap_height;
std::vector<Vec3f> envmap;


struct Light {
  Light(const Vec3f &p, const float i) : position(p), intensity(i) {}
  Vec3f position;
  float intensity;
};

struct Material {
  Material(const float r, const Vec4f &a, const Vec3f &color, const float spec) : refractive_index(r), albedo(a), diffuse_color(color), specular_exponent(spec) {}
  Material() : refractive_index(1), albedo(1,0,0,0), diffuse_color(), specular_exponent() {}
  float refractive_index;
  Vec4f albedo;
  Vec3f diffuse_color;
  float specular_exponent;
};



struct Sphere {
  Vec3f center;
  float radius;
  Material material;
  
  Sphere(const Vec3f &c, const float r, const Material &m) : center(c), radius(r), material(m) {}
  
  bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const {
    Vec3f L = center - orig;
    float tca = L*dir;
    float d2 = L*L - tca*tca;
    if (d2 > radius*radius) return false;
    float thc = sqrtf(radius*radius - d2);
    t0       = tca - thc;
    float t1 = tca + thc;
    if (t0 < 0) t0 = t1;
    if (t0 < 0) return false;
    return true;
  }
};

Vec3f reflect(const Vec3f &I, const Vec3f &N) {
  return I - N*2.f*(I*N);
}

Vec3f refract(const Vec3f &I, const Vec3f &N, const float eta_t, const float eta_i=1.f) { // Snell's law
  float cosi = - std::max(-1.f, std::min(1.f, I*N));
  if (cosi<0) return refract(I, -N, eta_i, eta_t); // if the ray comes from the inside the object, swap the air and the media
  float eta = eta_i / eta_t;
  float k = 1 - eta*eta*(1 - cosi*cosi);
  return k<0 ? Vec3f(1,0,0) : I*eta + N*(eta*cosi - sqrtf(k)); // k<0 = total reflection, no ray to refract. I refract it anyways, this has no physical meaning
}

bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, const std::vector<std::pair<Model, Material> > &models, Vec3f &hit, Vec3f &N, Material &material) {
  float dist = std::numeric_limits<float>::max();
  
  float spheres_dist = std::numeric_limits<float>::max();
  
  float checkerboard_dist = std::numeric_limits<float>::max();
  if (fabs(dir.y)>1e-3)  {
    float d = -(orig.y+4)/dir.y; // the checkerboard plane has equation y = -4
    Vec3f pt = orig + dir*d;
    if (d>0 && fabs(pt.x)<10 && pt.z<-10 && pt.z>-30 && d<spheres_dist) {
      checkerboard_dist = d;
      hit = pt;
      N = Vec3f(0,1,0);
      material.diffuse_color = (int(.5*hit.x+1000) + int(.5*hit.z)) & 1 ? Vec3f(.3, .3, .3) : Vec3f(.3, .2, .1);
    }
  }

  dist = std::min(dist, checkerboard_dist);

  for (size_t i = 0; i < models.size(); i++){ // itère sur tous les modèles
    Model model = models[i].first;
    for (int j = 0; j < model.nfaces(); j++){ 
      float dist_face;
      if ((model.ray_triangle_intersect(j, orig, dir, dist_face)) && (dist_face < dist)) {
        dist = dist_face;
        hit = orig + dir * dist_face;

        material = models[i].second; // on récupère le matériau du modèle
        
        Vec3f vert1 = model.point(model.vert(i,0));
        Vec3f vert2 = model.point(model.vert(i,1));
        Vec3f vert3 = model.point(model.vert(i,2));

        Vec3f edge1 = vert2 - vert1;
        Vec3f edge2 = vert3 - vert1;
        
        N = cross(edge1, edge2).normalize();
	
      }
    }
  }
  
  return dist < 1000;
}




Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, const std::vector<std::pair<Model, Material>> &models, const std::vector<Light> &lights, size_t depth=0) {
  Vec3f point, N;
  Material material;
  
  if (depth>4 || !scene_intersect(orig, dir, spheres, models, point, N, material)) {  
    int a = std::max(0, std::min(envmap_width -1, static_cast<int>((atan2(dir.z, dir.x)/(2*M_PI) + .5)*envmap_width)));
    int b = std::max(0, std::min(envmap_height-1, static_cast<int>(acos(dir.y)/M_PI*envmap_height)));
    return envmap[a+b*envmap_width];
  }
  
  Vec3f reflect_dir = reflect(dir, N).normalize();
  Vec3f refract_dir = refract(dir, N, material.refractive_index).normalize();
  Vec3f reflect_orig = reflect_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; // offset the original point to avoid occlusion by the object itself
  Vec3f refract_orig = refract_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
  Vec3f reflect_color = cast_ray(reflect_orig, reflect_dir, spheres, models, lights, depth + 1);
  Vec3f refract_color = cast_ray(refract_orig, refract_dir, spheres, models, lights, depth + 1);

  float diffuse_light_intensity = 0, specular_light_intensity = 0;
  for (size_t i=0; i<lights.size(); i++) {
    Vec3f light_dir      = (lights[i].position - point).normalize();
    float light_distance = (lights[i].position - point).norm();

    Vec3f shadow_orig = light_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; // checking if the point lies in the shadow of the lights[i]
    Vec3f shadow_pt, shadow_N;
    Material tmpmaterial;
    if (scene_intersect(shadow_orig, light_dir, spheres, models, shadow_pt, shadow_N, tmpmaterial) && (shadow_pt-shadow_orig).norm() < light_distance)
      continue;
	
    diffuse_light_intensity  += lights[i].intensity * std::max(0.f, light_dir*N);
    specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N)*dir), material.specular_exponent)*lights[i].intensity;
  }
  return material.diffuse_color * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.)*specular_light_intensity * material.albedo[1] + reflect_color*material.albedo[2] + refract_color*material.albedo[3];
}

void render(const std::vector<Sphere> &spheres, const std::vector<std::pair<Model, Material>> &models ,const std::vector<Light> &lights) {
  const int   width    = 1024;
  const int   height   = 768;
  const float fov      = M_PI/3.;
  std::vector<Vec3f> framebuffer(width*height);

#pragma omp parallel for
  for (size_t j = 0; j<height; j++) { // actual rendering loop
    for (size_t i = 0; i<width; i++) {
      float dir_x =  (i + 0.5) -  width/2.;
      float dir_y = -(j + 0.5) + height/2.;    // this flips the image at the same time
      float dir_z = -height/(2.*tan(fov/2.));

      Vec3f rendu;
      rendu = cast_ray(Vec3f(0,0,0), Vec3f(dir_x, dir_y, dir_z).normalize(), spheres, models, lights);
	    
      framebuffer[i+j*width] = (rendu);
    }
  }

  std::vector<unsigned char> pixmap(width*height*3);
  for (size_t i = 0; i < height*width; ++i) {
    Vec3f &c = framebuffer[i];
    float max = std::max(c[0], std::max(c[1], c[2]));
    if (max>1) c = c*(1./max);
    for (size_t j = 0; j<3; j++) {
      pixmap[i*3+j] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
    }
  }
  stbi_write_jpg("out.jpg", width, height, 3, pixmap.data(), 100);
}

int main() {
  int n = -1;
  unsigned char *pixmap = stbi_load("../snowy_forest.jpg", &envmap_width, &envmap_height, &n, 0);
  if (!pixmap || 3!=n) {
    std::cerr << "Error: can not load the environment map" << std::endl;
    return -1;
  }
  envmap = std::vector<Vec3f>(envmap_width*envmap_height);
  for (int j = envmap_height-1; j>=0 ; j--) {
    for (int i = 0; i<envmap_width; i++) {
      envmap[i+j*envmap_width] = Vec3f(pixmap[(i+j*envmap_width)*3+0], pixmap[(i+j*envmap_width)*3+1], pixmap[(i+j*envmap_width)*3+2])*(1/255.);
    }
  }
  stbi_image_free(pixmap);
  
  //Création des différents matériaux
  Material snow(1.0, Vec4f(0.4,  0.025, 0, 0.0), Vec3f(0.9, 0.9, 0.9),   10.);
  Material carrot_material(1.0, Vec4f(0.9,  0.1, 0.0, 0.0), Vec3f(0.9, 0.5, 0.1),   10.);
  Material wood(1.0, Vec4f(0.9,  0.02, 0.0, 0.0), Vec3f(0.25, 0.15, 0.08),   10.);
  Material button_eyes_material(1.0, Vec4f(0.9,  0.21, 0.0, 0.0), Vec3f(0.05, 0.05, 0.05),   10.);
  Material hat_soft(1.0, Vec4f(0.9,  0.025, 0.0, 0.0), Vec3f(0.11, 0.11, 0.11),   10.);

  std::vector<Sphere> spheres;
  std::vector<Light>  lights;
  lights.push_back(Light(Vec3f(-20, 20,  20), 1.5));
  lights.push_back(Light(Vec3f( 30, 20, -25), 1.8));
  lights.push_back(Light(Vec3f( 0, 20, 0), 1.7));  
  lights.push_back(Light(Vec3f( 30, 20, 30), 1.7));  
  lights.push_back(Light(Vec3f( 30, -20, 30), 1.7));  

  Model carrot = Model("../carrot.obj"); // on charge le modèle de la carotte
  
  std::pair<Model, Material> carrot_model = {carrot, carrot_material};
  
  std::vector<std::pair<Model, Material> > models;
  models.push_back(carrot_model);

  Model bodyb = Model("../bodyb.obj"); // on charge le modèle des boules de neiges
  std::pair<Model, Material> bodyb_model = {bodyb, snow};

  models.push_back(bodyb_model);

  Model hat = Model("../hat.obj"); // on charge le modèle du chapeau
  std::pair<Model, Material> hat_model = {hat, hat_soft};

  models.push_back(hat_model);

    Model eyes_button= Model("../eyes_button.obj"); //on charge le modèle des yeux et des boutons
  std::pair<Model, Material> eyes_buttons_model = {eyes_button, button_eyes_material};

  models.push_back(eyes_buttons_model);

  Model arms = Model("../arms.obj"); //on charge le modèle des bras
  std::pair<Model, Material> arms_model = {arms, wood};
  models.push_back(arms_model);

  render(spheres, models, lights);

  return 0;
}
