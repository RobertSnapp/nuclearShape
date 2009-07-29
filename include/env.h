#ifndef __env_h__
#define __env_h__
/******************************************************************************
 *
 *	File:			env.h
 *
 *	Function:		Select the opengl and glut headers that match the architecture.
 *
 *	Author:			Robert R. Snapp
 *
 *	Copyright:		Copyright (C) 2006 Robert R. Snapp
 *					All Rights Reserved.
 *
 *	Source:			Original.
 *
 *	Notes:			
 *
 *	Change History:
 *			11/16/06	Started source.
 *	
 ******************************************************************************/

#pragma once

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <Glut/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#endif /* __env_h__ */
