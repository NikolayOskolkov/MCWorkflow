#!/usr/bin/env python3
"""
Download files from Figshare.
Usage: python api.py <project_id> [target_file]
"""

import sys
import os
import requests
from pathlib import Path


def download_figshare_file(project_id, target_file=None):
    """
    Download a file from Figshare

    Args:
        project_id: Figshare project ID
        target_file: filename

    Returns:
        Path to downloaded file, or None if failed
    """
    # Figshare ndownloader API
    url = f"https://ndownloader.figshare.com/files/{project_id}"
    
    # Set up headers to mimic a browser and handle redirects
    headers = {
        'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36'
    }
    
    try:
        print(f"Connecting to {url}...")
        
        # Make request with redirects enabled and stream for progress
        response = requests.get(url, headers=headers, allow_redirects=True, stream=True, timeout=30)
        response.raise_for_status()
        
        # Get file size
        total_size = int(response.headers.get('content-length', 0))
        
        if total_size == 0:
            print("ERROR: Server returned 0 bytes. File may not exist or be inaccessible.")
            return None
        
        # Determine output filename
        if target_file is None:
            # Try to get filename from Content-Disposition header
            content_disp = response.headers.get('content-disposition', '')
            if 'filename=' in content_disp:
                target_file = content_disp.split('filename=')[-1].strip('"\'')
            else:
                target_file = f"figshare_file_{project_id}"

        output_path = Path(target_file)

        # Download with progress
        downloaded = 0
        print(f"Downloading to {target_file} ({total_size / (1024**3):.2f} GB)...")
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    # Simple progress indicator
                    if total_size > 0:
                        percent = (downloaded / total_size) * 100
                        if downloaded % (50 * 1024 * 1024) == 0:  # Print every 50MB
                            print(f"  {percent:.1f}% ({downloaded / (1024**3):.2f} GB)...")

        # Verify file was actually downloaded
        if output_path.stat().st_size == 0:
            print("ERROR: Downloaded file is empty!")
            output_path.unlink()
            return None

        print(f"✓ Download complete: {target_file} ({output_path.stat().st_size / (1024**3):.2f} GB)")
        return str(output_path)

    except requests.exceptions.RequestException as e:
        print(f"ERROR: Download failed: {e}")
        return None
    except Exception as e:
        print(f"ERROR: {e}")
        return None


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python api.py <project_id> [target_file]")
        print("Example: python api.py 52692896 GTDB_sliced_seqs_sliding_window.fna.gz")
        sys.exit(1)

    project_id = sys.argv[1]
    target_file = sys.argv[2] if len(sys.argv) > 2 else None

    result = download_figshare_file(project_id, target_file)
    sys.exit(0 if result else 1)
