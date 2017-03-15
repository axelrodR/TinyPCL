--
-- premake5 file to build RecastDemo
-- http://premake.github.io/
--

local action = _ACTION or ""
local todir = "Build/" .. action

solution "TinyPCL"
	configurations { 
		"Debug",
		"Release"
	}
	location (todir)

	-- extra warnings, no exceptions or rtti
	flags { 
		"ExtraWarnings",
		"FloatFast",
		"Symbols"
	}
	exceptionhandling "Off"
	rtti "Off"

	-- debug configs
	configuration "Debug*"
		defines { "DEBUG" }
		targetdir ( todir .. "/lib/Debug" )
 
 	-- release configs
	configuration "Release*"
		defines { "NDEBUG" }
		flags { "Optimize" }
		targetdir ( todir .. "/lib/Release" )

	-- windows specific
	configuration "windows"
		defines { "WIN32", "_WINDOWS", "_CRT_SECURE_NO_WARNINGS" }


project "TinyPCL"
	language "C++"
	kind "StaticLib"
	includedirs { 
		"../include",
		"../src",
		"../DetourTileCache/Include",
		"../Recast/Include"
	}
	files { 
		"*.h",
		"*.hpp",
		"*.cpp"
	}



project "Tests"
	language "C++"
	kind "ConsoleApp"

	-- Catch requires RTTI and exceptions
	exceptionhandling "On"
	rtti "On"

	includedirs { 
		"../test",
	}
	files	{ 
		"*.h",
		"*.hpp",
		"*.cpp",
	}

	-- project dependencies
	links { 
		"TinyPCL",
	}

	-- distribute executable in RecastDemo/Bin directory
	targetdir "bin"
