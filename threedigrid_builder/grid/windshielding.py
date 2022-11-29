import numpy as np

from threedigrid_builder.base import Array, search
from threedigrid_builder.constants import ContentType

__all__ = ["Windshieldings"]


class Windshielding:
    id: int
    channel_id: int
    north: float
    northeast: float
    east: float
    southeast: float
    south: float
    southwest: float
    west: float
    northwest: float


class Windshieldings(Array[Windshielding]):
    def apply_to_lines(self, lines):
        """Set windshielding factors on lines. Only channels can have windshielding,
        as only channels can have wind applied to them

            Args:
                lines (Lines)
        """
        if len(self) == 0:
            return

        line_ch_idx = np.where(lines.content_type == ContentType.TYPE_V2_CHANNEL)[0]
        has_windshielding = line_ch_idx[
            np.isin(lines.content_pk[line_ch_idx], self.channel_id)
        ]

        chn_idx = search(
            self.channel_id,
            lines.content_pk[has_windshielding],
            check_exists=False,
            assume_ordered=False,
        )

        lines.windshieldings[has_windshielding, 0] = self.north[chn_idx]
        lines.windshieldings[has_windshielding, 1] = self.northeast[chn_idx]
        lines.windshieldings[has_windshielding, 2] = self.east[chn_idx]
        lines.windshieldings[has_windshielding, 3] = self.southeast[chn_idx]
        lines.windshieldings[has_windshielding, 4] = self.south[chn_idx]
        lines.windshieldings[has_windshielding, 5] = self.southwest[chn_idx]
        lines.windshieldings[has_windshielding, 6] = self.west[chn_idx]
        lines.windshieldings[has_windshielding, 7] = self.northwest[chn_idx]
