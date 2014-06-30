Goals
-----
- Fill the void between computing SFHs (`match_wrapper`) and constructing a
  a final mosaic (`montage_wrapper`).

  - FSPS/`python-fsps`, `scombine`, and `sedpy` are used to process a given
    SFH. The gridding still needs to be done manually, hence the desire for a
    new package.
  - Abstracting the gridding process would allow images to be constructed from
    arbitrary data processed in an arbitrary way.

The new package should satisfy the following:

- Specify a list of files or a list/array of data in memory (each element in
  the list represents a cell in the grid; elements should be ordered as a
  flattened array from a FITS file).
- Specify a function that processes the data in a cell and produces a value for
  the output grid.

  - More generally, the function could take arbitrary arguments and keywords,
    i.e., ``func(*args, **kwargs)``. The input list could then be a list of
    ``args`` lists, and an additional list of ``kwargs`` could also be given.
    They would default to ``[[], [], ...]`` and ``[{}, {}, ...]``.

- By enforcing well-ordered input, a grid can be constructed simply by
  reshaping the output list of cell values. The number of rows and columns
  would have to be specified at minimum. Could also allow some options for
  flexibility, such as reversing the rows.
- Create a FITS HDU object (header optional) and write a FITS file for the
  grid.


Asides
------
- Use `montage_wrapper.mosaic` for mosaic production. Will still need to do
  preprocessing (e.g., masking, unit conversion, density conversion),
  postprocessing (e.g., density conversion, background/foreground subtraction),
  and file copying/moving/renaming myself. This is specific to the data and
  project, though, and doesn't belong in a generic package.
- `calc_sed` and `calc_mag` are useful abstractions pulling together the tools
  in `python_fsps`, `scombine`, and `sedpy` into a single interface. Kind of
  like "high level functions" such as `mosaic` in `montage_wrapper`. *Should
  these be in their own package?* It would be nice to somehow associate it with
  `scombine`/`sedpy`.

  - Name ideas: `fluxutil`, `scombine_wrapper`, `scombine_hlf`
  - What if Ben renames `scombine`?

- Should wcs stuff go into its own package? `wcsutil` (name taken?)
