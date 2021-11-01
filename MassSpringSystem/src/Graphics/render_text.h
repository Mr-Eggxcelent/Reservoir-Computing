///////////////////////////////////////////////////////////////////////////////////////////////
//Joey de Vries  https://twitter.com/JoeyDeVriez. 
//All code samples, unless explicitly stated otherwise, are licensed under the terms of the CC BY - 
//NC 4.0 license as published by Creative Commons, either version 4 of the License, or (at your option) 
//any later version.
//See https : learnopengl.com/About for more information.
//Code can be found at https://learnopengl.com/Getting-started/Shaders
//////////////////////////////////////////////////////////////////////////////////////////////
#pragma once
#include "shader.h"

#include <map>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <ft2build.h>
#include FT_FREETYPE_H


class RenderText
{

    unsigned int VAO, VBO;
    unsigned int WIDTH, HEIGHT;

    /// Holds all state information relevant to a character as loaded using FreeType
    struct Character {
        unsigned int TextureID; // ID handle of the glyph texture
        glm::ivec2   Size;      // Size of glyph
        glm::ivec2   Bearing;   // Offset from baseline to left/top of glyph
        unsigned int Advance;   // Horizontal offset to advance to next glyph
    };

    std::map<GLchar, Character> Characters;
    std::shared_ptr<Shader> _shader;
    std::shared_ptr<Shader> _bill_shader;
    glm::mat4 _projection;
    Camera& _camera;
    glm::mat4 _view;


public:
    RenderText(Camera& camera,unsigned int& Width,unsigned int& Height)
        :WIDTH(Width),HEIGHT(Height),_camera(camera)
    {

    }

    void init_text_renderer()
    {
        _shader = std::make_shared<Shader>("src/Shader/text.vs", "src/Shader/text.fs");

        // need to adjust this shader to allow it to been seen form all angles
        _bill_shader= std::make_shared<Shader>("src/Shader/billboard.vs", "src/Shader/billboard.fs");
        // FreeType
        // --------
        FT_Library ft;

    // All functions return a value different than 0 whenever an error occurred
       if (FT_Init_FreeType(&ft))
       {
           std::cout << "ERROR::FREETYPE: Could not init FreeType Library" << std::endl;
       }

      // find path to font
      std::string font_name = "src/Shader/arial.ttf";
      if (font_name.empty())
      {
          std::cout << "ERROR::FREETYPE: Failed to load font_name" << std::endl;
      }
      
      // load font as face
      FT_Face face;
      if (FT_New_Face(ft, font_name.c_str(), 0, &face)) {
          std::cout << "ERROR::FREETYPE: Failed to load font" << std::endl;
      }

      else {
          // set size to load glyphs as
          FT_Set_Pixel_Sizes(face, 0, 48);
      
          // disable byte-alignment restriction
          glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
      
          // load first 128 characters of ASCII set
          for (unsigned char c = 0; c < 128; c++)
          {
              // Load character glyph 
              if (FT_Load_Char(face, c, FT_LOAD_RENDER))
              {
                  std::cout << "ERROR::FREETYTPE: Failed to load Glyph" << std::endl;
                  continue;
              }
              // generate texture
              unsigned int texture;
              glGenTextures(1, &texture);
              glBindTexture(GL_TEXTURE_2D, texture);
              glTexImage2D(
                  GL_TEXTURE_2D,
                  0,
                  GL_RED,
                  face->glyph->bitmap.width,
                  face->glyph->bitmap.rows,
                  0,
                  GL_RED,
                  GL_UNSIGNED_BYTE,
                  face->glyph->bitmap.buffer
              );
              // set texture options
              glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
              glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
              glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
              glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
              // now store character for later use
              Character character = {
                  texture,
                  glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
                  glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
                  static_cast<unsigned int>(face->glyph->advance.x)
              };
              Characters.insert(std::pair<char, Character>(c, character));
          }
          glBindTexture(GL_TEXTURE_2D, 0);
      }
      // destroy FreeType once we're finished
      FT_Done_Face(face);
      FT_Done_FreeType(ft);
      
      
      // configure VAO/VBO for texture quads
      // -----------------------------------
      glGenVertexArrays(1, &VAO);
      glGenBuffers(1, &VBO);
      glBindVertexArray(VAO);
      glBindBuffer(GL_ARRAY_BUFFER, VBO);
      glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 6 * 4, NULL, GL_DYNAMIC_DRAW);
      glEnableVertexAttribArray(0);
      glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), 0);
      glBindBuffer(GL_ARRAY_BUFFER, 0);
      glBindVertexArray(0);
      
    }
//    // render line of text
//// -------------------
    void render_billboard(std::string text, float x, float y, float scale, glm::vec3 color)
    {    
        glm::mat4 projection = glm::perspective(glm::radians(_camera.Zoom), (float)WIDTH / (float)HEIGHT, 0.1f, 100.0f);
        _view = _camera.GetViewMatrix();

        _bill_shader->use();
        glUniformMatrix4fv(glGetUniformLocation(_bill_shader->GetID(), "projection"), 1, GL_FALSE, glm::value_ptr(projection));
        glUniformMatrix4fv(glGetUniformLocation(_bill_shader->GetID(), "view"), 1, GL_FALSE, glm::value_ptr(_view));
        _bill_shader->setVec3("CameraRight_worldspace", glm::vec3(_camera.GetViewMatrix()[0][0], _camera.GetViewMatrix()[1][0], _camera.GetViewMatrix()[2][0]));
        _bill_shader->setVec3("CameraUp_worldspace", glm::vec3(_camera.GetViewMatrix()[0][1], _camera.GetViewMatrix()[1][1], _camera.GetViewMatrix()[2][1]));

        glUniform3f(glGetUniformLocation(_bill_shader->GetID(), "textColor"), color.x, color.y, color.z);
        glActiveTexture(GL_TEXTURE0);
        glBindVertexArray(VAO);

        // iterate through all characters
        std::string::const_iterator c;
        for (c = text.begin(); c != text.end(); c++)
        {
            Character ch = Characters[*c];

            float xpos = x + ch.Bearing.x * scale;
            float ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

            float w = ch.Size.x * scale;
            float h = ch.Size.y * scale;
            // update VBO for each character
            float vertices[6][4] = {
                { xpos,     ypos + h,   0.0f, 0.0f },
                { xpos,     ypos,       0.0f, 1.0f },
                { xpos + w, ypos,       1.0f, 1.0f },

                { xpos,     ypos + h,   0.0f, 0.0f },
                { xpos + w, ypos,       1.0f, 1.0f },
                { xpos + w, ypos + h,   1.0f, 0.0f }
            };

            _bill_shader->setVec3("center_pos", glm::vec3(xpos,ypos,0));
           
            // render glyph texture over quad
            glBindTexture(GL_TEXTURE_2D, ch.TextureID);
            // update content of VBO memory
            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); // be sure to use glBufferSubData and not glBufferData

            glBindBuffer(GL_ARRAY_BUFFER, 0);
            // render quad
            glDrawArrays(GL_TRIANGLES, 0, 6);
            // now advance cursors for next glyph (note that advance is number of 1/64 pixels)
            x += (ch.Advance >> 6) * scale; // bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
        }
        glBindVertexArray(0);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    void render_text(std::string text, float x, float y, float scale, glm::vec3 color)
    {
        _projection = glm::ortho(0.0f, static_cast<float>(WIDTH), 0.0f, static_cast<float>(HEIGHT));
        _shader->use();
        glUniformMatrix4fv(glGetUniformLocation(_shader->GetID(), "projection"), 1, GL_FALSE, glm::value_ptr(_projection));
        glUniform3f(glGetUniformLocation(_shader->GetID(), "textColor"), color.x, color.y, color.z);
        glActiveTexture(GL_TEXTURE0);
        glBindVertexArray(VAO);

        // iterate through all characters
        std::string::const_iterator c;
        for (c = text.begin(); c != text.end(); c++)
        {
            Character ch = Characters[*c];

            float xpos = x + ch.Bearing.x * scale;
            float ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

            float w = ch.Size.x * scale;
            float h = ch.Size.y * scale;
            // update VBO for each character
            float vertices[6][4] = {
                { xpos,     ypos + h,   0.0f, 0.0f },
                { xpos,     ypos,       0.0f, 1.0f },
                { xpos + w, ypos,       1.0f, 1.0f },

                { xpos,     ypos + h,   0.0f, 0.0f },
                { xpos + w, ypos,       1.0f, 1.0f },
                { xpos + w, ypos + h,   1.0f, 0.0f }
            };
            // render glyph texture over quad
            glBindTexture(GL_TEXTURE_2D, ch.TextureID);
            // update content of VBO memory
            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); // be sure to use glBufferSubData and not glBufferData

            glBindBuffer(GL_ARRAY_BUFFER, 0);
            // render quad
            glDrawArrays(GL_TRIANGLES, 0, 6);
            // now advance cursors for next glyph (note that advance is number of 1/64 pixels)
            x += (ch.Advance >> 6) * scale; // bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
        }
        glBindVertexArray(0);
        glBindTexture(GL_TEXTURE_2D, 0);
    }



};