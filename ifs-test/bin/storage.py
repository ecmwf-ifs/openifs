#! /usr/bin/env python3

import logging
import pathlib
import re
import shutil
import urllib.error
import urllib.request

import click
import yaml

class DataSource:
  """
  Base class for data sources.
  """

  def available(self):
    """
    Check whether the given data source is accessible.
    """
    return NotImplemented

  def fetch(self, source_path, target_path):
    """
    Fetch a file.

    :param source_path:     The relative path to the source file.
    :param target_path:     The absolute path where the targe file should be
                            placed.
    """
    return NotImplemented


class FileSource(DataSource):
  """
  Use a given (cache) path on the local filesystem as the data source.
  """

  def __init__(self, root):
    self._root_path = pathlib.Path(root).resolve()

  def available(self):
    return self._root_path.exists()

  def fetch(self, source_path, target_path):
    source_path = self._root_path/source_path

    if not source_path.exists():
      raise RuntimeError("The source file %s does not exist!" % source_path)

    target_path.parent.mkdir(parents=True, exist_ok=True)
    target_path.unlink(missing_ok=True)

    logging.debug(f"Create symlink from {source_path} to {target_path}")
    target_path.symlink_to(source_path)

  def __str__(self):
    return "FileSource(%s)" % self._root_path


class DownloadSource(DataSource):
  """
  Download a file from a given URL.
  """

  def __init__(self, root):
    self._url = str(root)


  def available(self):
    try:
      urllib.request.urlopen(self._url)
    except urllib.error.URLError:
      return False

    return True

  def fetch(self, source_path, target_path):
    target_path.parent.mkdir(parents=True, exist_ok=True)
    target_path.unlink(missing_ok=True)

    url = self._url + "/" + source_path
    logging.debug(f"Download file from {url} to {target_path}.")

    try:
      urllib.request.urlretrieve(url, filename=target_path)
    except urllib.error.URLError as ue:
      raise RuntimeError(str(ue))

  def __str__(self):
    return "DownloadSource(%s)" % self._url

def _open_config(path):
  """
  Read the configuration from a yaml file and return it as a dictionary.
  """
  with open(path, 'r') as f:
    ref_data = yaml.safe_load(f)

  return ref_data

def _get_data_sources(config):
  """
  From the configuration, find the available data sources.

  :param config:    The configuration as a dictionary.
  """

  sources = []

  for data in config['external_data']['data_source']:
    kwargs = dict(data)
    kwargs.pop('type')

    source = None

    if data['type'] == 'file':
      source = FileSource(**kwargs)
    elif data['type'] == 'download':
      source = DownloadSource(**kwargs)

    if source is not None:
      if source.available():
        sources.append(source)
        logging.debug(f"Source {source} is available.")
      else:
        logging.debug(f"Source {source} is not available.")

  return sources

@click.group()
@click.option('-v', 'verbosity', count=True)
def cli(verbosity):
  if verbosity == 0:
    level = logging.CRITICAL
  elif verbosity == 1:
    level = logging.INFO
  else:
    level = logging.DEBUG

  logging.basicConfig(level=level)

@cli.command()
@click.argument('config_file',
  type=click.Path(exists=True, path_type=pathlib.Path, dir_okay=False),
)
@click.argument('output_dir',
  type=click.Path(path_type=pathlib.Path, file_okay=False),
)
@click.option('--force/--no-force',
  default=False,
  help="Overwrite existing files"
)
@click.option('-p', '--custom-path',
  type=click.Path(file_okay=False, dir_okay=True, exists=True),
  multiple=True,
  help="Add custom data source path."
)
def fetch(config_file, output_dir, force, custom_path):
  """
  Fetch all files that are specified in the given config YAML file.

  Each specified file will be placed in `output_dir/target_path` where
  `target_path` is the corresponding target location that is specified in the
  config file.

  If `extract_path` is also given in the config file, the archive will be extracted
  to `output_dir/extract_path`.

  Parameters
  ----------
  config_file:
    The configuration YAML file that holds all file specifications.
  output_dir:
    The directory where all the fetched files are put.
  """
  logging.info("Start fetching data.")
  logging.debug("Configuration file = %s." % config_file)
  logging.debug("Output directory = %s." % output_dir)
  logging.debug("Force = %s." % force)

  config = _open_config(config_file)
  sources = _get_data_sources(config)

  if custom_path:
    sources = [FileSource(path) for path in custom_path] + sources

  if not sources:
    raise RuntimeError("No accessible data source was found.")

  for data in config['external_data']['files']:
    source_path = data['source_path']
    target_path = data['target_path']
    extract_path = data.get('extract_path', None)

    # target_path is a relative path. To get the absolute path, it must be
    # combined with the target (cache) directory.
    target_path = output_dir/target_path

    if target_path.exists() and force:
      logging.debug("Delete %s." % target_path)
      target_path.unlink()

    if target_path.exists():
      logging.debug(f"Skip file {target_path} (file exists).")
    else:
      for source in sources:
        success = False
        try:
          source.fetch(source_path, target_path)
          success = True
          break
        except RuntimeError as re:
          logging.debug(f"Fetching {source_path} from {source} failed.")

      if not success:
        raise RuntimeError(f"Fetching {source_path} failed!")

    if extract_path is not None:
      extract_path = output_dir/extract_path

      if extract_path.exists() and force:
        logging.debug("Delete %s." % extract_path)
        shutil.rmtree(extract_path)

      if extract_path.exists():
        logging.debug(f"Skip extracting to {extract_path} (directory exists).")
      else:
        shutil.unpack_archive(target_path, extract_path)

  logging.info("Fetching data finished. Data has been placed in %s." % output_dir)

@cli.command()
@click.argument('config_file',
  type=click.Path(exists=True, path_type=pathlib.Path, dir_okay=False),
)
@click.argument('output_dir',
  type=click.Path(path_type=pathlib.Path, file_okay=False),
)
@click.option('--force/--no-force',
  default=False,
  help="Overwrite existing files"
)
@click.option('-p', '--custom-path',
  type=click.Path(file_okay=False, dir_okay=True, exists=True),
  multiple=True,
  help="Add custom data source path."
)
def cache(config_file, output_dir, force, custom_path):
  """
  Fetch all files that are specified in the given config YAML file and put
  them in `output_dir`, using the same directory structure as in the data 
  source.

  This means that each specified file will be placed in `output_dir/source_path`
  where `source_path` is the corresponding target location that is specified in
  the config file.

  Parameters
  ----------
  config_file:
    The configuration YAML file that holds all file specifications.
  output_dir:
    The directory where all the fetched files are put.

  """
  logging.info("Start caching data.")
  logging.debug("Configuration file = %s." % config_file)
  logging.debug("Output directory = %s." % output_dir)
  logging.debug("Force = %s." % force)

  config = _open_config(config_file)
  sources = _get_data_sources(config)

  if custom_path:
    sources = [FileSource(path) for path in custom_path] + sources

  if not sources:
    raise RuntimeError("No accessible data source was found.")

  for data in config['external_data']['files']:
    source_path = data['source_path']

    target_path = output_dir/source_path

    if target_path.exists() and force:
      logging.debug("Delete %s." % target_path)
      target_path.unlink()

    if target_path.exists():
      logging.debug(f"Skip file {target_path} (file exists).")
    else:
      for source in sources:
        success = False
        try:
          source.fetch(source_path, target_path)
          success = True
          break
        except RuntimeError as re:
          logging.debug(f"Fetching {source_path} from {source} failed.")

      if not success:
        raise RuntimeError(f"Fetching {source_path} failed!")

  logging.info("Caching data finished. Data has been placed in %s." % output_dir)

@cli.command()
@click.argument('config_file',
  type=click.Path(exists=True, path_type=pathlib.Path, dir_okay=False),
)
@click.argument('cache_dir',
  type=click.Path(path_type=pathlib.Path, file_okay=False),
)
@click.argument('output_dir',
  type=click.Path(path_type=pathlib.Path, file_okay=False),
)
@click.option('--force/--no-force',
  default=False,
  help="Overwrite existing files"
)
@click.option('--match',
  default="",
  help="Only process entries where source_path matches this Python-style regex"
)
def pack(config_file, cache_dir, output_dir, force, match):
  """
  Repack archives that are specified in the config YAML file.

  For each file that has `extract_path` specified in the config file, create an
  archive from `cache_dir/extract_path` and put it into
  `output_dir/source_path`.

  Parameters
  ----------
  config_file:
    The configuration YAML file that holds all file specifications.
  cache_dir:
    The root directory where the extracted archives are located. Usually 
    corresponds to `output_dir` in `storage.py fetch`.
  output_dir:
    The directory where the archives will be put.
  """
  logging.info("Start packing data.")
  logging.debug("Configuration file = %s." % config_file)
  logging.debug("Cache directory = %s." % cache_dir)
  logging.debug("Output directory = %s." % output_dir)
  logging.debug("Force = %s." % force)

  config = _open_config(config_file)

  for data in config['external_data']['files']:
    source_path = data['source_path']
    target_path = data['target_path']
    extract_path = data.get('extract_path', None)

    # Ignore all non-archive entries (no extract_path specified).
    if extract_path is None:
      continue

    if match and (re.search(match, source_path) is None):
      logging.info(f"Skip {source_path} as it doesn't match {match}.")
      continue

    # Path to the directory which we want to archive.
    data_dir = cache_dir/extract_path

    # Path of the archive that we want to create.
    archive_path = output_dir/source_path

    if not data_dir.exists():
      logging.info(f"Directory {data_dir} does not exist. Won't build {target_path}.")

    if archive_path.exists():
      if force:
        logging.debug("Delete {archive_path}.")
        archive_path.unlink()
      else:
        logging.info(f"File {archive_path} exists already. Skipping...")
        continue

    # Get the archive suffix without the .
    archive_type = archive_path.suffix[1:]
    if archive_type == 'gz':
      archive_type = 'gztar'

    logging.debug(f"Pack {data_dir} into {archive_path}")

    # Ugly suffix hacking. As make_archive always adds the archive type suffix,
    # we have to remove it first. We do .stem.stem to remove double-suffixes
    # like tar.gz.
    shutil.make_archive(
      archive_path.parent/pathlib.Path(archive_path.stem).stem,
      archive_type,
      data_dir
    )

  logging.info("Packing data finished.")


if __name__ == "__main__":
  cli()
