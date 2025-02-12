def arg_to_file_exist(arg, mandatory=True, io="in", err = 0):

  import os

  filename = None
  if arg == None:
    if mandatory:
      print("ERROR: mandatory argument missing")
      err = 1
  else:
      filename=arg
      exists = os.path.isfile(filename)

      if (exists) & (io == "out"):
        print("ERROR: output file already exists ",filename)
        err = 2
      if (not exists) & (io == "in"):
        if mandatory:
          print("ERROR: input file does not exist ",filename)
          err = 3
        else:
          print("WARNING: input file does not exist ",filename)
          filename = None

  return filename, err
