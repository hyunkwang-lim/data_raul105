import os
import ycm_core

database = ycm_core.CompilationDatabase('build')

def FlagsForFile(filename, **kwargs):
  if database:
    path_without_ext = os.path.splitext(filename)[0]
    candidate_files = [
        filename,
        path_without_ext + '.c',
        path_without_ext + '.cpp',
        os.path.abspath('src/input/drivers/segmenter/main/main.cpp'),
        ]
    for candidate_file in candidate_files:
      compilation_info = database.GetCompilationInfoForFile(candidate_file)
      if compilation_info.compiler_flags_: break

    if compilation_info:
      return {
        'flags': ['-xc++'] + list(compilation_info.compiler_flags_),
        'do_cache': True
      }

  return { 'flags': [] }
