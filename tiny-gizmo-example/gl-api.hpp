// This is free and unencumbered software released into the public domain.
// For more information, please refer to <http://unlicense.org>

#pragma once

#ifndef gl_api_hpp
#define gl_api_hpp

#include <map>
#include <string>
#include <vector>
#include "linalg.h"

namespace
{
	static bool gEnableGLDebugOutputErrorBreakpoints = false;

	inline void compile_shader(GLuint program, GLenum type, const char * source)
	{
		GLuint shader = glCreateShader(type);
		glShaderSource(shader, 1, &source, nullptr);
		glCompileShader(shader);

		GLint status, length;
		glGetShaderiv(shader, GL_COMPILE_STATUS, &status);

		if (status == GL_FALSE)
		{
			glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);
			std::vector<GLchar> buffer(length);
			glGetShaderInfoLog(shader, (GLsizei)buffer.size(), nullptr, buffer.data());
			glDeleteShader(shader);
			std::cerr << "GL Compile Error: " << buffer.data() << std::endl;
			std::cerr << "Source: " << source << std::endl;
			throw std::runtime_error("GLSL Compile Failure");
		}

		glAttachShader(program, shader);
		glDeleteShader(shader);
	}

	std::string gl_src_to_str(GLenum source)
	{
		switch (source)
		{
		case GL_DEBUG_SOURCE_WINDOW_SYSTEM: return "WINDOW_SYSTEM";
		case GL_DEBUG_SOURCE_SHADER_COMPILER: return "SHADER_COMPILER";
		case GL_DEBUG_SOURCE_THIRD_PARTY: return "THIRD_PARTY";
		case GL_DEBUG_SOURCE_APPLICATION: return "APPLICATION";
		case GL_DEBUG_SOURCE_OTHER: return "OTHER";
		case GL_DEBUG_SOURCE_API: return "API";
		default: return "UNKNOWN";
		}
	}

	std::string gl_enum_to_str(GLenum type)
	{
		switch (type)
		{
		case GL_DEBUG_TYPE_ERROR: return "ERROR";
		case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR: return "DEPRECATION";
		case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR: return "UNDEFINED_BEHAVIOR";
		case GL_DEBUG_TYPE_PORTABILITY: return "PORTABILITY";
		case GL_DEBUG_TYPE_PERFORMANCE: return "PERFORMANCE";
		case GL_DEBUG_TYPE_OTHER: return "OTHER";
		default: return "UNKNOWN";
		}
	}

	std::string gl_severity_to_str(GLenum severity)
	{
		switch (severity)
		{
		case GL_DEBUG_SEVERITY_LOW: return "LOW";
		case GL_DEBUG_SEVERITY_MEDIUM: return "MEDIUM";
		case GL_DEBUG_SEVERITY_HIGH: return "HIGH";
		default: return "UNKNOWN";
		}
	}

	static void APIENTRY gl_debug_callback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar * message, GLvoid * userParam)
	{
		if (type != GL_DEBUG_TYPE_ERROR) return;
		auto sourceStr = gl_src_to_str(source);
		auto typeStr = gl_enum_to_str(type);
		auto severityStr = gl_severity_to_str(severity);
		std::cout << "gl_debug_callback: " << sourceStr << ", " << severityStr << ", " << typeStr << " , " << id << ", " << message << std::endl;
		if ((type == GL_DEBUG_TYPE_ERROR) && (gEnableGLDebugOutputErrorBreakpoints)) __debugbreak();
	}

	inline void gl_check_error(const char * file, int32_t line)
	{
#if defined(_DEBUG) || defined(DEBUG)
		GLint error = glGetError();
		if (error)
		{
			const char * errorStr = 0;
			switch (error)
			{
			case GL_INVALID_ENUM: errorStr = "GL_INVALID_ENUM"; break;
			case GL_INVALID_VALUE: errorStr = "GL_INVALID_VALUE"; break;
			case GL_INVALID_OPERATION: errorStr = "GL_INVALID_OPERATION"; break;
			case GL_OUT_OF_MEMORY: errorStr = "GL_OUT_OF_MEMORY"; break;
			default: errorStr = "unknown error"; break;
			}
			printf("GL error : %s, line %d : %s\n", file, line, errorStr);
			error = 0;
		}
#endif
	}

	inline size_t gl_size_bytes(GLenum type)
	{
		switch (type)
		{
		case GL_UNSIGNED_BYTE: return sizeof(uint8_t);
		case GL_UNSIGNED_SHORT: return sizeof(uint16_t);
		case GL_UNSIGNED_INT: return sizeof(uint32_t);
		default: throw std::logic_error("unknown element type"); break;
		}
	}
}

template<typename factory_t>
class GlObject
{
	mutable GLuint handle = 0;
	std::string n;
public:
	GlObject() {}
	GlObject(GLuint h) : handle(g) {}
	~GlObject() { if (handle) factory_t::destroy(handle); }
	GlObject(const GlObject & r) = delete;
	GlObject & operator = (GlObject && r) { std::swap(handle, r.handle); std::swap(n, r.n);  return *this; }
	GlObject(GlObject && r) { *this = std::move(r); }
	operator GLuint () const { if (!handle)  factory_t::create(handle); return handle; }
	GlObject & operator = (GLuint & other) { handle = other; return *this; }
	void set_name(const std::string & newName) { n = newName; }
	std::string name() const { return n; }
	GLuint id() const { return handle; };
};

struct GlBufferFactory { static void create(GLuint & x) { glGenBuffers(1, &x); }; static void destroy(GLuint x) { glDeleteBuffers(1, &x); }; };
struct GlTextureFactory { static void create(GLuint & x) { glGenTextures(1, &x); }; static void destroy(GLuint x) { glDeleteTextures(1, &x); }; };
struct GlVertexArrayFactory { static void create(GLuint & x) { glGenVertexArrays(1, &x); }; static void destroy(GLuint x) { glDeleteVertexArrays(1, &x); }; };
struct GlRenderbufferFactory { static void create(GLuint & x) { glGenRenderbuffers(1, &x); }; static void destroy(GLuint x) { glDeleteRenderbuffers(1, &x); }; };
struct GlFramebufferFactory { static void create(GLuint & x) { glGenFramebuffers(1, &x); }; static void destroy(GLuint x) { glDeleteFramebuffers(1, &x); }; };
struct GlQueryFactory { static void create(GLuint & x) { glGenQueries(1, &x); }; static void destroy(GLuint x) { glDeleteQueries(1, &x); }; };
struct GlSamplerFactory { static void create(GLuint & x) { glGenSamplers(1, &x); }; static void destroy(GLuint x) { glDeleteSamplers(1, &x); }; };
struct GlTransformFeedbacksFactory { static void create(GLuint & x) { glGenTransformFeedbacks(1, &x); }; static void destroy(GLuint x) { glDeleteTransformFeedbacks(1, &x); }; };

typedef GlObject<GlBufferFactory> GlBufferObject;
typedef GlObject<GlTextureFactory> GlTextureObject;
typedef GlObject<GlVertexArrayFactory> GlVertexArrayObject;
typedef GlObject<GlRenderbufferFactory> GlRenderbufferObject;
typedef GlObject<GlFramebufferFactory> GlFramebufferObject;
typedef GlObject<GlQueryFactory> GlQueryObject;
typedef GlObject<GlSamplerFactory> GlSamplerObject;
typedef GlObject<GlTransformFeedbacksFactory> GlTransformFeedbacksObject;

//////////////////
//   GlBuffer   //
//////////////////

struct GlBuffer : public GlBufferObject
{
	GLsizeiptr size;
	GlBuffer() {}
	void set_buffer_data(const GLsizeiptr s, const GLvoid * data, const GLenum usage) { this->size = s; glNamedBufferDataEXT(*this, size, data, usage);  }
	void set_buffer_data(const std::vector<GLubyte> & bytes, const GLenum usage) { set_buffer_data(bytes.size(), bytes.data(), usage); }
    void set_buffer_sub_data(const GLsizeiptr s, const GLintptr offset, const GLvoid * data) { glNamedBufferSubDataEXT(*this, offset, s, data);  }
    void set_buffer_sub_data(const std::vector<GLubyte> & bytes, const GLintptr offset, const GLenum usage) { set_buffer_sub_data(bytes.size(), offset, bytes.data()); }
};

////////////////////////
//   GlRenderbuffer   //
////////////////////////

struct GlRenderbuffer : public GlRenderbufferObject
{
	float width{ 0 }, height{ 0 };
	GlRenderbuffer() {}
	GlRenderbuffer(float width, float height) : width(width), height(height) {}
};

///////////////////////
//   GlFramebuffer   //
///////////////////////

struct GlFramebuffer : public GlFramebufferObject
{
	float width{ 0 }, height{ 0 }, depth{ 0 };
	GlFramebuffer() {}
	GlFramebuffer(float width, float height) : width(width), height(height) {}
	GlFramebuffer(float width, float height, float depth) : width(width), height(height), depth(depth) {}
	void check_complete() { if (glCheckNamedFramebufferStatusEXT(*this, GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) throw std::runtime_error("fbo incomplete"); }
};

///////////////////
//   GlTexture   //
///////////////////

struct GlTexture2D : public GlTextureObject
{
	float width{ 0 }, height{ 0 };
	GlTexture2D() {}
	GlTexture2D(float width, float height) : width(width), height(height) {}

	void setup(GLsizei width, GLsizei height, GLenum internal_fmt, GLenum format, GLenum type, const GLvoid * pixels, bool createMipmap = false)
	{
		glTextureImage2DEXT(*this, GL_TEXTURE_2D, 0, internal_fmt, width, height, 0, format, type, pixels);
		if (createMipmap) glGenerateTextureMipmapEXT(*this, GL_TEXTURE_2D);
		glTextureParameteriEXT(*this, GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTextureParameteriEXT(*this, GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, createMipmap ? GL_LINEAR_MIPMAP_LINEAR : GL_LINEAR);
		glTextureParameteriEXT(*this, GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTextureParameteriEXT(*this, GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		this->width = width;
		this->height = height;
	}
};

/////////////////////
//   GlTexture3D   //
/////////////////////

// As either a 3D texture or 2D array
struct GlTexture3D : public GlTextureObject
{
	float width{ 0 }, height{ 0 }, depth{ 0 };
    GlTexture3D() {}
	GlTexture3D(float width, float height, float depth) : width(width), height(height), depth(depth) {}

	void setup(GLenum target, GLsizei width, GLsizei height, GLsizei depth, GLenum internal_fmt, GLenum format, GLenum type, const GLvoid * pixels)
	{
		glTextureImage3DEXT(*this, target, 0, internal_fmt, width, height, depth, 0, format, type, pixels);
		glTextureParameteriEXT(*this, target, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTextureParameteriEXT(*this, target, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTextureParameteriEXT(*this, target, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
		glTextureParameteriEXT(*this, target, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
		glTextureParameteriEXT(*this, target, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER);
		this->width = width;
		this->height = height;
		this->depth = depth;
	}
};

//////////////////
//   GlShader   //
//////////////////

class GlShader
{
	GLuint program;
	bool enabled = false;

protected:
	GlShader(const GlShader & r) = delete;
	GlShader & operator = (const GlShader & r) = delete;
public:

	GlShader() : program() {}

    GlShader(const GLuint type, const std::string & src)
    {
        program = glCreateProgram();

        ::compile_shader(program, type, src.c_str());
        glProgramParameteri(program, GL_PROGRAM_SEPARABLE, GL_TRUE);

        glLinkProgram(program);

        GLint status, length;
        glGetProgramiv(program, GL_LINK_STATUS, &status);

        if (status == GL_FALSE)
        {
            glGetProgramiv(program, GL_INFO_LOG_LENGTH, &length);
            std::vector<GLchar> buffer(length);
            glGetProgramInfoLog(program, (GLsizei)buffer.size(), nullptr, buffer.data());
            std::cerr << "GL Link Error: " << buffer.data() << std::endl;
            throw std::runtime_error("GLSL Link Failure");
        }
    }

	GlShader(const std::string & vert, const std::string & frag, const std::string & geom = "")
	{
		program = glCreateProgram();

        glProgramParameteri(program, GL_PROGRAM_SEPARABLE, GL_FALSE);

		::compile_shader(program, GL_VERTEX_SHADER, vert.c_str());
		::compile_shader(program, GL_FRAGMENT_SHADER, frag.c_str());

		if (geom.length() != 0) ::compile_shader(program, GL_GEOMETRY_SHADER, geom.c_str());

		glLinkProgram(program);

		GLint status, length;
		glGetProgramiv(program, GL_LINK_STATUS, &status);

		if (status == GL_FALSE)
		{
			glGetProgramiv(program, GL_INFO_LOG_LENGTH, &length);
			std::vector<GLchar> buffer(length);
			glGetProgramInfoLog(program, (GLsizei)buffer.size(), nullptr, buffer.data());
			std::cerr << "GL Link Error: " << buffer.data() << std::endl;
			throw std::runtime_error("GLSL Link Failure");
		}
	}

	~GlShader() { if (program) glDeleteProgram(program); }

	GlShader(GlShader && r) : GlShader() { *this = std::move(r); }

	GLuint handle() const { return program; }
	GLint get_uniform_location(const std::string & name) const { return glGetUniformLocation(program, name.c_str()); }

	GlShader & operator = (GlShader && r) { std::swap(program, r.program); return *this; }

	std::map<uint32_t, std::string> reflect()
	{
        std::map<uint32_t, std::string> locations;
		GLint count;
		glGetProgramiv(program, GL_ACTIVE_UNIFORMS, &count);
		for (GLuint i = 0; i < static_cast<GLuint>(count); ++i)
		{
			char buffer[1024]; GLenum type; GLsizei length; GLint size, block_index;
			glGetActiveUniform(program, i, sizeof(buffer), &length, &size, &type, buffer);
			glGetActiveUniformsiv(program, 1, &i, GL_UNIFORM_BLOCK_INDEX, &block_index);
			if (block_index != -1) continue;
			GLint loc = glGetUniformLocation(program, buffer);
            locations[loc] = std::string(buffer);
		}
        return locations;
	}

	void uniform(const std::string & name, int scalar) const { glProgramUniform1i(program, get_uniform_location(name), scalar); }
	void uniform(const std::string & name, float scalar) const { glProgramUniform1f(program, get_uniform_location(name), scalar); }
	void uniform(const std::string & name, const linalg::aliases::float2 & vec) const { glProgramUniform2fv(program, get_uniform_location(name), 1, &vec.x); }
	void uniform(const std::string & name, const linalg::aliases::float3 & vec) const { glProgramUniform3fv(program, get_uniform_location(name), 1, &vec.x); }
	void uniform(const std::string & name, const linalg::aliases::float4 & vec) const { glProgramUniform4fv(program, get_uniform_location(name), 1, &vec.x); }
	void uniform(const std::string & name, const linalg::aliases::float3x3 & mat) const { glProgramUniformMatrix3fv(program, get_uniform_location(name), 1, GL_FALSE, &mat.x.x); }
	void uniform(const std::string & name, const linalg::aliases::float4x4 & mat) const { glProgramUniformMatrix4fv(program, get_uniform_location(name), 1, GL_FALSE, &mat.x.x); }

	void uniform(const std::string & name, const int elements, const std::vector<int> & scalar) const { glProgramUniform1iv(program, get_uniform_location(name), elements, scalar.data()); }
	void uniform(const std::string & name, const int elements, const std::vector<float> & scalar) const { glProgramUniform1fv(program, get_uniform_location(name), elements, scalar.data()); }
	void uniform(const std::string & name, const int elements, const std::vector<linalg::aliases::float2> & vec) const { glProgramUniform2fv(program, get_uniform_location(name), elements, &vec[0].x); }
	void uniform(const std::string & name, const int elements, const std::vector<linalg::aliases::float3> & vec) const { glProgramUniform3fv(program, get_uniform_location(name), elements, &vec[0].x); }
	void uniform(const std::string & name, const int elements, const std::vector<linalg::aliases::float3x3> & mat) const { glProgramUniformMatrix3fv(program, get_uniform_location(name), elements, GL_FALSE, &mat[0].x.x); }
	void uniform(const std::string & name, const int elements, const std::vector<linalg::aliases::float4x4> & mat) const { glProgramUniformMatrix4fv(program, get_uniform_location(name), elements, GL_FALSE, &mat[0].x.x); }

	void texture(GLint loc, GLenum target, int unit, GLuint tex) const
	{
		glBindMultiTextureEXT(GL_TEXTURE0 + unit, target, tex);
		glProgramUniform1i(program, loc, unit);
	}

	void texture(const char * name, int unit, GLuint tex, GLenum target) const { texture(get_uniform_location(name), target, unit, tex); }

	void bind() { if (program > 0) enabled = true; glUseProgram(program); }
	void unbind() { enabled = false; glUseProgram(0); }
};

////////////////
//   GlMesh   //
////////////////

class GlMesh
{
	GlVertexArrayObject vao;
	GlBuffer vertexBuffer, instanceBuffer, indexBuffer;

	GLenum drawMode = GL_TRIANGLES;
	GLenum indexType = 0;
	GLsizei vertexStride = 0, instanceStride = 0, indexCount = 0;

public:

	GlMesh() {}
	GlMesh(GlMesh && r) { *this = std::move(r); }
	GlMesh(const GlMesh & r) = delete;
	GlMesh & operator = (GlMesh && r)
	{
		char buffer[sizeof(GlMesh)];
		memcpy(buffer, this, sizeof(buffer));
		memcpy(this, &r, sizeof(buffer));
		memcpy(&r, buffer, sizeof(buffer));
		return *this;
	}
	GlMesh & operator = (const GlMesh & r) = delete;
	~GlMesh() {};

	void set_non_indexed(GLenum newMode)
	{
		drawMode = newMode;
		indexBuffer = {};
		indexType = 0;
		indexCount = 0;
	}

	void draw_elements(int instances = 0) const
	{
		if (vertexBuffer.size)
		{
			glBindVertexArray(vao);
			if (indexCount)
			{
				glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
				if (instances) glDrawElementsInstanced(drawMode, indexCount, indexType, 0, instances);
				else glDrawElements(drawMode, indexCount, indexType, nullptr);
			}
			else
			{
				if (instances) glDrawArraysInstanced(drawMode, 0, static_cast<GLsizei>(vertexBuffer.size / vertexStride), instances);
				else glDrawArrays(drawMode, 0, static_cast<GLsizei>(vertexBuffer.size / vertexStride));
			}
			glBindVertexArray(0);
		}
	}

	void set_vertex_data(GLsizeiptr size, const GLvoid * data, GLenum usage) { vertexBuffer.set_buffer_data(size, data, usage); }
    GlBuffer & get_vertex_data_buffer() { return vertexBuffer; };

	void set_instance_data(GLsizeiptr size, const GLvoid * data, GLenum usage) { instanceBuffer.set_buffer_data(size, data, usage); }

	void set_index_data(GLenum mode, GLenum type, GLsizei count, const GLvoid * data, GLenum usage)
	{
		size_t size = gl_size_bytes(type);
		indexBuffer.set_buffer_data(size * count, data, usage);
		drawMode = mode;
		indexType = type;
		indexCount = count;
	}
    GlBuffer & get_index_data_buffer() { return indexBuffer; };

	void set_attribute(GLuint index, GLint size, GLenum type, GLboolean normalized, GLsizei stride, const GLvoid * offset)
	{
		glEnableVertexArrayAttribEXT(vao, index);
		glVertexArrayVertexAttribOffsetEXT(vao, vertexBuffer, index, size, type, normalized, stride, (GLintptr)offset);
		vertexStride = stride;
	}

	void set_instance_attribute(GLuint index, GLint size, GLenum type, GLboolean normalized, GLsizei stride, const GLvoid * offset)
	{
		glEnableVertexArrayAttribEXT(vao, index);
		glVertexArrayVertexAttribOffsetEXT(vao, instanceBuffer, index, size, type, normalized, stride, (GLintptr)offset);
		glVertexArrayVertexAttribDivisorEXT(vao, index, 1);
		instanceStride = stride;
	}

	void set_indices(GLenum mode, GLsizei count, const uint8_t * indices, GLenum usage) { set_index_data(mode, GL_UNSIGNED_BYTE, count, indices, usage); }
	void set_indices(GLenum mode, GLsizei count, const uint16_t * indices, GLenum usage) { set_index_data(mode, GL_UNSIGNED_SHORT, count, indices, usage); }
	void set_indices(GLenum mode, GLsizei count, const uint32_t * indices, GLenum usage) { set_index_data(mode, GL_UNSIGNED_INT, count, indices, usage); }

	template<class T> void set_vertices(size_t count, const T * vertices, GLenum usage) { set_vertex_data(count * sizeof(T), vertices, usage); }
	template<class T> void set_vertices(const std::vector<T> & vertices, GLenum usage) { set_vertices(vertices.size(), vertices.data(), usage); }
	template<class T, int N> void set_vertices(const T(&vertices)[N], GLenum usage) { set_vertices(N, vertices, usage); }

	template<class V> void set_attribute(GLuint index, float V::*field) { set_attribute(index, 1, GL_FLOAT, GL_FALSE, sizeof(V), &(((V*)0)->*field)); }
	template<class V, int N> void set_attribute(GLuint index, linalg::vec<float, N> V::*field) { set_attribute(index, N, GL_FLOAT, GL_FALSE, sizeof(V), &(((V*)0)->*field)); }

	template<class T> void set_elements(GLsizei count, const linalg::vec<T, 2> * elements, GLenum usage) { set_indices(GL_LINES, count * 2, &elements->x, usage); }
	template<class T> void set_elements(GLsizei count, const linalg::vec<T, 3> * elements, GLenum usage) { set_indices(GL_TRIANGLES, count * 3, &elements->x, usage); }
	template<class T> void set_elements(GLsizei count, const linalg::vec<T, 4> * elements, GLenum usage) { set_indices(GL_QUADS, count * 4, &elements->x, usage); }

	template<class T> void set_elements(const std::vector<T> & elements, GLenum usage) { set_elements((GLsizei)elements.size(), elements.data(), usage); }

	template<class T, int N> void set_elements(const T(&elements)[N], GLenum usage) { set_elements(N, elements, usage); }
};

#endif // end gl_api_hpp
